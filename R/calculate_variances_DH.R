#' Segregation variance for DH lines from specified crosses
#'
#' Calculate expected genomic estimated breeding values (EGBV), segregation
#' variance, and superior progeny value (SPV) for doubled haploid (DH) lines
#' derived after \code{t} rounds of random mating.
#'
#' @param crosses A data.frame or matrix (n_crosses x 2) of parental indices
#'   (row indices into \code{M}) specifying the crosses.
#' @param M Numeric marker matrix (markers in columns, genotypes in rows) with
#'   entries in \code{c(0, 2)} corresponding to \code{c("AA", "BB")}.
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(M))
#'     \item \code{chr}: chromosome identifier (numeric)
#'     \item \code{pos}: position on the chromosome in morgan
#'   }
#'   All markers in \code{M} must appear in \code{genetic.map$site}.
#' @param U Numeric matrix of marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(U) == ncol(M)}.
#' @param t Integer. Number of random‐mating generations before DH creation.
#' @param intensity Double. Standardized selection differential.
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param calculate.index Logical. If \code{TRUE}, calculate the index from fixed
#'   trait \code{weights}. Only works if \code{covariance = TRUE}; otherwise ignored.
#' @param weights Numeric vector of length \code{ncol(U)} with fixed trait weights.
#' @param method Character. Which method to use; one of
#'   \code{c("osthushenrich", "lehermeier")}.
#' @param n.Threads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named \code{EGBV1, var1, SPV1, EGBV2, ...} and a list with element \code{covariances}.
#'   If \code{covariance = TRUE} and \code{calculate.index = TRUE}, additionally an element \code{index}
#'   (a data.frame) with \code{IDG_A, VARIDG_A, SPVIDG_A}.
#'   If \code{covariance = FALSE}, a data.frame with the same columns \code{EGBV1, var1, SPV1, EGBV2, ...}.
#' @export
calculate_variances_DH <- function(crosses, genetic.map, M, U, t, intensity,
                                   covariance = FALSE, calculate.index = FALSE, weights = NULL,
                                   method = "osthushenrich", n.Threads = 4L) {

  # ---- Normalize ----
  if (!is.matrix(M)) M <- as.matrix(M)
  if (!is.matrix(U)) U <- as.matrix(U)
  crosses <- as.matrix(crosses)
  if(ncol(U)==1 & is.null(weights)){weights <- c(1)}
  # Validate method early
  method <- match.arg(method, c("osthushenrich", "lehermeier"))

  # ---- Checks ----
  if (ncol(M) <= 0L) stop("M must have markers in columns.")
  if (nrow(U) != ncol(M)) stop("U must have nrow(U) == ncol(M).")
  if (!is.logical(covariance) || length(covariance) != 1L) stop("`covariance` must be TRUE/FALSE.")
  if (!is.numeric(t) || length(t) != 1L || t < 0 || abs(t - round(t)) > .Machine$double.eps^0.5)
    stop("`t` must be a single non-negative integer.")
  if (!is.numeric(intensity) || length(intensity) != 1L)
    stop("`intensity` must be a single numeric.")
  if (ncol(crosses) != 2L) stop("`crosses` must have 2 columns.")

  # Weights checks (only required if index requested)
  if (calculate.index && is.null(weights)) {
    stop("`weights` is required when calculate.index = TRUE.")
  }
  if (calculate.index && length(weights) != ncol(U)) {
    stop("`weights` must have length equal to ncol(U) when calculate.index = TRUE.")
  }
  if (calculate.index) {
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) stop("`weights` must contain only finite values.")
  }

  # t must be >= 1 (you used t<0 before; that was a bug)
  if (t <= 0) {
    message("t needs to be >= 1, setting t = 1")
    t <- 1
  }

  # Validate crosses
  genos <- unique(as.vector(crosses))
  if (any(!is.finite(genos))) stop("`crosses` contains non-finite entries.")
  if (any(genos < 1 | genos > nrow(M))) stop("Some genotype indices in `crosses` are outside 1..nrow(M).")

  # Validate map
  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map)))
    stop("`genetic.map` must contain columns: site, chr, pos.")
  if (any(genetic.map$site < 1 | genetic.map$site > ncol(M)))
    stop("`genetic.map$site` contains indices outside 1..ncol(M).")
  if (!all(seq_len(ncol(M)) %in% genetic.map$site))
    stop("Some markers in M are missing from `genetic.map$site`.")

  # ---- Reorder M and U to (chr, pos) like in your 4W wrapper (this was missing before) ----
  map2 <- genetic.map[order(genetic.map$chr, genetic.map$pos), , drop = FALSE]
  ord  <- as.integer(map2$site)

  M <- M[, ord, drop = FALSE]
  U <- U[ord,  , drop = FALSE]

  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  # One-trait shortcut
  if (ncol(U) == 1L) covariance <- FALSE

  # If covariance disabled, index not available
  if (!covariance && calculate.index) {
    message("covariance == FALSE: ignoring index related arguments")
    calculate.index <- FALSE
    weights <- rep(0, ncol(U))
  }

  # Threads
  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  # ---- Call C++ ----
  # IMPORTANT: this assumes your C++ functions now accept:
  #   weights + calcindex (not gains + calcgains)
  if (method == "lehermeier") {
    temp <- cpp_calculate_covariance_lehermeier(
      Crosses    = crosses,
      genMap     = genmap_list,
      M          = M,
      U          = U,
      t          = as.integer(t),
      intensity  = intensity,
      weights    = weights,
      covariance = covariance,
      calcindex  = calculate.index,
      nThreads   = nThreads
    )
  } else {
    temp <- cpp_calculate_covariance_osthushenrich(
      Crosses    = crosses,
      genMap     = genmap_list,
      M          = M,
      U          = U,
      t          = as.integer(t),
      intensity  = intensity,
      weights    = weights,
      covariance = covariance,
      calcindex  = calculate.index,
      nThreads   = nThreads
    )
  }

  # ---- Format outputs ----
  name_vec <- paste0(rep(c("EGBV","var","SPV"), each = ncol(U)), seq_len(ncol(U)))
  crosses_df <- as.data.frame(crosses)
  names(crosses_df) <- c("parent1", "parent2")

  if (covariance) {
    if (calculate.index) {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(3 * ncol(U)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)

      temp2 <- as.data.frame(temp$cross_values)[, ((3 * ncol(U)) + 1):ncol(as.data.frame(temp$cross_values)), drop = FALSE]
      names(temp2) <- c("IDG_A", "VARIDG_A", "SPVIDG_A")
      temp2 <- cbind(crosses_df, temp2)

      return(list(cross_values = temp1, index = temp2, covariances = temp$covariances))
    } else {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(3 * ncol(U)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)
      return(list(cross_values = temp1, covariances = temp$covariances))
    }
  } else {
    out <- as.data.frame(temp)
    names(out) <- name_vec
    out <- cbind(crosses_df, out)
    return(out)
  }
}
