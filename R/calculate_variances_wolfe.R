#' Additive and dominance segregation variance of families from specified crosses (Wolfe)
#'
#' Calculate expected genomic estimated breeding value (EGBV), expected total genetic value (ETGV),
#' additive and dominance segregation variance, superior progeny value with additive variance (SPV),
#' and superior progeny value including dominance variance (TSPV) for F1 families.
#'
#' @param crosses A data.frame or matrix (n_crosses x 2) of parental indices
#'   (row indices into \code{hap.mat1}/\code{hap.mat2}) specifying the crosses.
#' @param hap.mat1 Numeric haplotype matrix (individuals x markers) for the first haplotype,
#'   entries in \code{c(0, 1)} for \code{c("A", "B")}.
#' @param hap.mat2 Numeric haplotype matrix (individuals x markers) for the second haplotype,
#'   entries in \code{c(0, 1)} for \code{c("A", "B")}.
#' @param genetic.map A data.frame with columns:
#'   \itemize{
#'     \item \code{site}: integer marker index (1..ncol(hap.mat1))
#'     \item \code{chr}: chromosome identifier (numeric or factor)
#'     \item \code{pos}: position on the chromosome in morgan
#'   }
#'   All markers in \code{hap.mat1}/\code{hap.mat2} must appear in \code{genetic.map$site}.
#' @param effects.A Numeric matrix of additive marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(effects.A) == ncol(hap.mat1)}.
#' @param effects.D Numeric matrix of dominance marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(effects.D) == ncol(hap.mat1)} and \code{ncol(effects.D) == ncol(effects.A)}.
#' @param intensity Double. Standardized selection differential (used for SPV).
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param weights Numeric vector of length \code{ncol(effects.A)} with fixed trait weights.
#' If supplied the function calculates index values for each cross. Only works if \code{covariance = TRUE}; otherwise ignored.
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
calculate_variances_F1 <- function(crosses, genetic.map, hap.mat1, hap.mat2, effects.A, effects.D,
                                   intensity, covariance, weights = NULL,
                                   n.Threads = 4L) {

  hap.mat1 <- as.matrix(hap.mat1); hap.mat2 <- as.matrix(hap.mat2)
  effects.A <- as.matrix(effects.A); effects.D <- as.matrix(effects.D)
  crosses <- as.matrix(crosses)
  calculate.index <- !is.null(weights)

  if (ncol(effects.A) == 1 & is.null(weights)) { weights <- c(1) }

  if (!is.logical(covariance) || length(covariance) != 1L) {
    stop("`covariance` must be a single logical (TRUE/FALSE).")
  }
  if (!is.numeric(intensity) || length(intensity) != 1L) {
    stop("`intensity` must be a single numeric (standardized selection differential).")
  }

  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  if (ncol(hap.mat1) <= 0L) stop("hap.mat1 must have markers in columns.")
  if (!identical(dim(hap.mat1), dim(hap.mat2))) stop("hap.mat1 and hap.mat2 must have identical dimensions.")
  if (nrow(effects.A) != ncol(hap.mat1)) {
    stop("effects.A must have nrow(effects.A) == ncol(hap.mat1). Found: nrow(effects.A) = ", nrow(effects.A), ", ncol(hap.mat1) = ", ncol(hap.mat1), ".")
  }
  if (nrow(effects.D) != ncol(hap.mat1)) {
    stop("effects.D must have nrow(effects.D) == ncol(hap.mat1). Found: nrow(effects.D) = ", nrow(effects.D), ", ncol(hap.mat1) = ", ncol(hap.mat1), ".")
  }
  if (ncol(effects.D) != ncol(effects.A)) stop("effects.A and effects.D must have the same number of trait columns.")
  if (ncol(crosses) != 2L) stop("`crosses` must have exactly 2 columns (P1, P2).")

  if (calculate.index && is.null(weights)) {
    stop("`weights` is required when calculate.index = TRUE.")
  }
  if (calculate.index && length(weights) != ncol(effects.A)) {
    stop("`weights` must have length equal to ncol(effects.A) when calculate.index = TRUE.")
  }
  if (calculate.index) {
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) stop("`weights` must contain only finite values.")
  }

  if (!covariance && calculate.index) {
    message("covariance == FALSE: ignoring index related arguments")
    calculate.index <- FALSE
    weights <- rep(0, ncol(effects.A))
  }

  genos <- unique(as.vector(crosses))
  if (any(!is.finite(genos))) stop("`crosses` contains non-finite entries.")
  if (any(genos < 1 | genos > nrow(hap.mat1))) {
    stop("Some genotype indices in `crosses` are outside 1..nrow(hap.mat1).")
  }

  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map))) {
    stop("`genetic.map` must contain columns: ", paste(req_cols, collapse = ", "), ".")
  }
  if (any(genetic.map$site < 1 | genetic.map$site > ncol(hap.mat1))) {
    stop("`genetic.map$site` contains indices outside 1..ncol(hap.mat1).")
  }
  if (!all(seq_len(ncol(hap.mat1)) %in% genetic.map$site)) {
    stop("Some markers in hap.mat1/hap.mat2 are missing from `genetic.map$site`.")
  }

  map2 <- genetic.map[order(genetic.map$chr, genetic.map$pos), , drop = FALSE]
  ord  <- as.integer(map2$site)

  hap.mat1 <- hap.mat1[, ord, drop = FALSE]
  hap.mat2 <- hap.mat2[, ord, drop = FALSE]
  effects.A <- effects.A[ord, , drop = FALSE]
  effects.D <- effects.D[ord, , drop = FALSE]

  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  temp <- cpp_calculate_covariance_wolfe(
    Crosses    = crosses,
    genMap     = genmap_list,
    Hap1       = hap.mat1,
    Hap2       = hap.mat2,
    U          = effects.A,
    D          = effects.D,
    intensity  = intensity,
    weights    = weights,
    covariance = covariance,
    calcindex  = calculate.index,
    nThreads   = nThreads
  )

  name_vec <- paste0(rep(c("EGBV","ETGV","var_a","SPV","var_d","TSPV"), each = ncol(effects.A)), seq_len(ncol(effects.A)))

  crosses_df <- as.data.frame(crosses)
  names(crosses_df) <- c("parent1", "parent2")

  if (covariance) {
    if (calculate.index) {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(6 * ncol(effects.A)), drop = FALSE]
      names(temp1) <- name_vec
      temp1 <- cbind(crosses_df, temp1)

      temp2 <- as.data.frame(temp$cross_values)[, (6 * ncol(effects.A) + 1):(6 * ncol(effects.A) + 6), drop = FALSE]
      names(temp2) <- c("IDG_A","VARIDG_A","SPVIDG_A","IDG_AD","VARIDG_AD","SPVIDG_AD")
      temp2 <- cbind(crosses_df, temp2)
      if (temp$check_psd) {
        warning("Some segregation covariance matrices were not psd; for index calculation the nearest psd was used in these places")
      }
      return(list(
        cross_values = temp1,
        index = temp2,
        additive_covariances = temp$covA,
        dominance_covariances = temp$covD
      ))
    } else {
      temp1 <- as.data.frame(temp$cross_values)[, 1:(6 * ncol(effects.A)), drop = FALSE]
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
