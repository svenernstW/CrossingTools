#' Additive and dominance segregation variance of families from specified crosses (Wolfe)
#'
#' Calculate expected genomic estimated breeding value (EGBV), expected total genetic value (ETGV),
#' additive and dominance segregation variance, superior progeny value with additive variance (SPV),
#' and superior progeny value including dominance variance (TSPV) for F1 families.
#'
#' @param crosses A data.frame or matrix (n_crosses x 2) of parental indices
#'   (row indices into \code{hap1}/\code{hap2}) specifying the crosses.
#' @param hap1 Numeric haplotype matrix (individuals x markers) for the first haplotype,
#'   entries in \code{c(0, 1)} for \code{c("A", "B")}.
#' @param hap2 Numeric haplotype matrix (individuals x markers) for the second haplotype,
#'   entries in \code{c(0, 1)} for \code{c("A", "B")}.
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(hap1))
#'     \item \code{chr}: chromosome identifier (numeric or factor)
#'     \item \code{pos}: position on the chromosome in morgan
#'   }
#'   All markers in \code{hap1}/\code{hap2} must appear in \code{genetic.map$site}.
#' @param U Numeric matrix of additive marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(U) == ncol(hap1)}.
#' @param D Numeric matrix of dominance marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(D) == ncol(hap1)} and \code{ncol(D) == ncol(U)}.
#' @param intensity Double. Standardized selection differential (used for SPV).
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param calculate.index Logical. If \code{TRUE}, calculate the index from fixed
#'   trait \code{weights}. Only works if \code{covariance = TRUE}; otherwise ignored.
#' @param weights Numeric vector of length \code{ncol(U)} with fixed trait weights.
#'
#' @param n.Threads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross_values}
#'   (a data.frame) whose columns are named
#'   \code{EGBV1, ETGV1, varA1, SPV1, varD1, TSPV1, EGBV2, ...}
#'   and covariance matrices per cross.
#'   If \code{covariance = TRUE} and \code{calculate.index = TRUE}, additionally an element \code{index}
#'   (a data.frame) with \code{IDG_A, VARIDG_A, SPVIDG_A, IDG_AD, VARIDG_AD, SPVIDG_AD}.
#'   If \code{covariance = FALSE}, a data.frame with the same columns.
#'
#' @export
calculate_variances_F1 <- function(crosses, genetic.map, hap1, hap2, U, D,
                                   intensity, covariance,
                                   calculate.index = FALSE, weights = NULL,
                                   n.Threads = 4L) {

  hap1 <- as.matrix(hap1); hap2 <- as.matrix(hap2)
  U <- as.matrix(U); D <- as.matrix(D)
  crosses <- as.matrix(crosses)
  if(ncol(U)==1 & is.null(weights)){weights <- c(1)}
  # ---- Basic checks ----
  if (!is.logical(covariance) || length(covariance) != 1L) {
    stop("`covariance` must be a single logical (TRUE/FALSE).")
  }
  if (!is.numeric(intensity) || length(intensity) != 1L) {
    stop("`intensity` must be a single numeric (standardized selection differential).")
  }

  # Threads
  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  if (ncol(hap1) <= 0L) stop("hap1 must have markers in columns.")
  if (!identical(dim(hap1), dim(hap2))) stop("hap1 and hap2 must have identical dimensions.")
  if (nrow(U) != ncol(hap1)) {
    stop("U must have nrow(U) == ncol(hap1). Found: nrow(U) = ", nrow(U), ", ncol(hap1) = ", ncol(hap1), ".")
  }
  if (nrow(D) != ncol(hap1)) {
    stop("D must have nrow(D) == ncol(hap1). Found: nrow(D) = ", nrow(D), ", ncol(hap1) = ", ncol(hap1), ".")
  }
  if (ncol(D) != ncol(U)) stop("U and D must have the same number of trait columns.")
  if (ncol(crosses) != 2L) stop("`crosses` must have exactly 2 columns (P1, P2).")

  # ---- Weights checks (strict, like your 4-way wrapper) ----
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

  # If covariance is FALSE, index is not available
  if (!covariance && calculate.index) {
    message("covariance == FALSE: ignoring index related arguments")
    calculate.index <- FALSE
    weights <- rep(0, ncol(U))
  }

  # ---- Validate crosses against haplotypes ----
  genos <- unique(as.vector(crosses))
  if (any(!is.finite(genos))) stop("`crosses` contains non-finite entries.")
  if (any(genos < 1 | genos > nrow(hap1))) {
    stop("Some genotype indices in `crosses` are outside 1..nrow(hap1).")
  }

  # ---- Validate and align genetic.map with markers ----
  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map))) {
    stop("`genetic.map` must contain columns: ", paste(req_cols, collapse = ", "), ".")
  }
  if (any(genetic.map$site < 1 | genetic.map$site > ncol(hap1))) {
    stop("`genetic.map$site` contains indices outside 1..ncol(hap1).")
  }
  if (!all(seq_len(ncol(hap1)) %in% genetic.map$site)) {
    stop("Some markers in hap1/hap2 are missing from `genetic.map$site`.")
  }

  # ---- Crucial: reorder to (chr, pos), and rebuild genmap consistently ----
  map2 <- genetic.map[order(genetic.map$chr, genetic.map$pos), , drop = FALSE]
  ord  <- as.integer(map2$site)

  hap1 <- hap1[, ord, drop = FALSE]
  hap2 <- hap2[, ord, drop = FALSE]
  U    <- U[ord, , drop = FALSE]
  D    <- D[ord, , drop = FALSE]

  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  # ---- Call C++ (must match your updated signature) ----
  temp <- cpp_calculate_covariance_wolfe(
    Crosses    = crosses,
    genMap     = genmap_list,
    Hap1       = hap1,
    Hap2       = hap2,
    U          = U,
    D          = D,
    intensity  = intensity,
    weights    = weights,
    covariance = covariance,
    calcindex  = calculate.index,
    nThreads   = nThreads
  )

  # ---- Format outputs ----
  name_vec <- paste0(rep(c("EGBV","ETGV","var_a","SPV","var_d","TSPV"), each = ncol(U)), seq_len(ncol(U)))

  crosses_df <- as.data.frame(crosses)
  names(crosses_df) <- c("parent1", "parent2")

  if (covariance) {
    if (calculate.index) {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(6 * ncol(U)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)

      temp2 <- as.data.frame(temp$cross_values)[, (6 * ncol(U) + 1):(6 * ncol(U) + 6), drop = FALSE]
      names(temp2) <- c("IDG_A","VARIDG_A","SPVIDG_A","IDG_AD","VARIDG_AD","SPVIDG_AD")
      temp2 <- cbind(crosses_df, temp2)

      return(list(
        cross_values = temp1,
        index = temp2,
        additive_covariances = temp$covA,
        dominance_covariances = temp$covD
      ))
    } else {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(6 * ncol(U)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)

      return(list(
        cross_values = temp1,
        additive_covariances = temp$covA,
        dominance_covariances = temp$covD
      ))
    }
  } else {
    out <- as.data.frame(temp)
    names(out) <- name_vec
    out <- cbind(crosses_df, out)
    return(out)
  }
}
