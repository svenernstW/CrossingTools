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
#' @param marker.effects.A Numeric matrix of additive marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(marker.effects.A) == ncol(hap.mat1)}.
#' @param marker.effects.D Numeric matrix of dominance marker effects (markers in rows, traits in columns).
#'   Must have \code{nrow(effects.D) == ncol(hap.mat1)} and \code{ncol(marker.effects.D) == ncol(marker.effects.A)}.
#' @param intensity Double. Standardized selection differential (used for SPV), defaults to 1.
#' @param weights Numeric vector of length \code{ncol(marker.effects.A)} with fixed trait weights.
#' If supplied the function calculates index values for each cross.
#' @param covariance Logical. If \code{TRUE}, also compute the segregation
#'   covariance matrix between all supplied traits.
#' @param nthreads Integer (default 4). Number of OpenMP threads (if enabled at compile time).
#'
#' @return If \code{covariance = TRUE}, a list with element \code{cross.df}
#'   (a data.frame) whose columns are named
#'   \code{GEBV1, TGV1, var.A1, SPV1, var.D1, TSPV1, EGEBV2, ...}
#'   and covariance matrices per cross.
#'   If weights are supplied, additionally an element \code{index.df}
#'   (a data.frame) with \code{IDX.A, var.IDX.A, SPV.IDX, IDX.AD, var.IDX.AD, SPV.IDX.AD}.
#'   If no weights are supplied, a data.frame with the same columns.
#'
#' @export
calc_spv_outcross <- function(crosses, genetic.map, hap.mat1, hap.mat2, marker.effects.A, marker.effects.D,
                                   intensity=NULL, weights = NULL, covariance = FALSE,
                                   nthreads = 4L) {
  if(!ncol(crosses) %in% c(2)){stop("ncol(crosses) needs to be 2 ")}
  crosses_in <- crosses

  hap.mat1 <- as.matrix(hap.mat1); hap.mat2 <- as.matrix(hap.mat2)
  effects.A <- as.matrix(marker.effects.A); effects.D <- as.matrix(marker.effects.D)
  crosses <- as.matrix(crosses)

  n.Threads <- nthreads

  if (is.numeric(crosses) || is.integer(crosses)) {
    if (any(!is.finite(crosses))) stop("`crosses` contains non-finite entries.")
    if (any(crosses < 1 | crosses > nrow(hap.mat1))) {
      stop("Some genotype indices in `crosses` are outside 1..nrow(hap.mat1).")
    }
    crosses2 <- as.data.frame(crosses)
  } else {
    if (is.null(rownames(hap.mat1))) stop("Character `crosses` requires rownames(hap.mat1).")

    idx <- match(as.vector(crosses), rownames(hap.mat1))
    if (anyNA(idx)) stop("Some entries in `crosses` are not in rownames(hap.mat1).")

    idx <- matrix(idx, nrow = nrow(crosses), ncol = ncol(crosses), byrow = FALSE)
    crosses2 <- as.data.frame(idx)
  }

  crosses2 <- as.matrix(crosses2)
  storage.mode(crosses2) <- "integer"
  if (anyNA(crosses2)) stop("Internal error: `crosses2` contains NA after conversion.")


  if (!is.logical(covariance) || length(covariance) != 1L) {
    stop("`covariance` must be a single logical (TRUE/FALSE).")
  }
  if (!is.numeric(intensity) || length(intensity) != 1L) {
    warning("`intensity` must be a single numeric (standardized selection differential), setting it to 1.")
    intensity <- 1
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

  ntraits <- ncol(effects.A)

  if (ntraits == 1L && !is.null(weights)) {
    warning("Single trait: weights are ignored (index equals trait).", call. = FALSE)
    weights <- NULL
  }
  if (ntraits == 1L) {
    covariance <- FALSE
  }

  calculate.index <- !is.null(weights) && covariance
  calculate.simple.index <- !is.null(weights) && !covariance

  if (is.null(weights)) {
    weights <- rep(1, ntraits)
  } else {
    weights <- as.numeric(weights)
    if (length(weights) != ntraits) {
      stop("`weights` must have length equal to ncol(effects.A).", call. = FALSE)
    }
    if (any(!is.finite(weights))) {
      stop("`weights` must contain only finite values.", call. = FALSE)
    }
  }

  if (calculate.simple.index) {
    effects.A <- cbind(effects.A, effects.A %*% weights)
    effects.D <- cbind(effects.D, effects.D %*% weights)
  }

  req_cols <- c("site", "chr", "pos")
  if (!all(req_cols %in% names(genetic.map))){
    stop("`genetic.map` must contain columns: site, chr, pos.")
  }
  site <- genetic.map$site

  # Only needed if site is character
  if (is.numeric(site) || is.integer(site)) {

    site_idx <- as.integer(site)

    if (any(site_idx < 1 | site_idx > ncol(hap.mat1), na.rm = TRUE))
      stop("`genetic.map$site` contains indices outside 1..ncol(hap.mat1).")

    if (!all(seq_len(ncol(hap.mat1)) %in% site_idx))
      stop("Some markers in hap.mat1 and 2 are missing from `genetic.map$site`.")

  } else {

    mnames <- colnames(hap.mat1)
    if (is.null(mnames))
      stop("Character `genetic.map$site` requires colnames(hap.mat1).")

    site_chr <- as.character(site)

    if (anyDuplicated(site_chr))
      stop("`genetic.map$site` contains duplicated marker names.")

    if (!all(mnames %in% site_chr))
      stop("Some markers in hap.mat1 are missing from `genetic.map$site`.")

    site_idx <- match(site_chr, mnames)
    if (anyNA(site_idx))
      stop("`genetic.map$site` contains names not found in colnames(hap.mat1).")
  }

  # Order map by chr/pos and carry indices along
  o    <- order(genetic.map$chr, genetic.map$pos)
  map2 <- genetic.map[o, , drop = FALSE]
  ord  <- site_idx[o]


  # relabel sites to 1..p for downstream code
  map2$site <- seq_len(ncol(hap.mat1))


  hap.mat1 <- hap.mat1[, ord, drop = FALSE]
  hap.mat2 <- hap.mat2[, ord, drop = FALSE]
  effects.A <- effects.A[ord, , drop = FALSE]
  effects.D <- effects.D[ord, , drop = FALSE]

  chr_levels <- unique(map2$chr)
  genmap_list <- lapply(chr_levels, function(cc) {
    as.matrix(map2$pos[map2$chr == cc])
  })

  temp <- cpp_calculate_covariance_wolfe(
    Crosses    = crosses2,
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


  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- c("parent1", "parent2")
  name_vec <- c(
    paste0("GEBV", seq_len(ntraits)),
    paste0("TGV", seq_len(ntraits)),
    paste0("var.A", seq_len(ntraits)),
    paste0("SPV", seq_len(ntraits)),
    paste0("var.D", seq_len(ntraits)),
    paste0("TSPV", seq_len(ntraits))
  )

  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- c("parent1", "parent2")

  if (covariance) {
    cv <- as.data.frame(temp$cross_values)

    temp1 <- cv[, 1:(6 * ntraits), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    if (calculate.index) {
      temp2 <- cv[, (6 * ntraits + 1):(6 * ntraits + 6), drop = FALSE]
      names(temp2) <- c("GEBV.IDX","TGV.IDX","var.A.IDX","SPV.IDX","var.D.IDX","TSPV.IDX")
      temp2 <- cbind(crosses_df, temp2)

      if (isTRUE(temp$check_psd)) {
        warning("Some segregation covariance matrices were not psd; for index calculation the nearest psd was used in these places.", call. = FALSE)
      }

      return(list(
        cross.df = temp1,
        index.df = temp2,
        additive.covariances = temp$covA,
        dominance.covariances = temp$covD
      ))
    } else {
      return(list(
        cross.df = temp1,
        additive.covariances = temp$covA,
        dominance.covariances = temp$covD
      ))
    }
  }

  if (calculate.simple.index) {
    cv <- as.data.frame(temp)

    p <- ncol(effects.A)  # original traits + 1 appended index trait

    trait_cols <- c(
      seq_len(ntraits),
      p + seq_len(ntraits),
      2 * p + seq_len(ntraits),
      3 * p + seq_len(ntraits),
      4 * p + seq_len(ntraits),
      5 * p + seq_len(ntraits)
    )

    idx_cols <- c(p, 2 * p, 3 * p, 4 * p, 5 * p, 6 * p)

    temp1 <- cv[, trait_cols, drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    temp2 <- cv[, idx_cols, drop = FALSE]
    names(temp2) <- c("GEBV.IDX","TGV.IDX","var.A.IDX","SPV.IDX","var.D.IDX","TSPV.IDX")
    temp2 <- cbind(crosses_df, temp2)

    return(list(
      cross.df = temp1,
      index.df = temp2
    ))
  }

  out <- as.data.frame(temp)
  out <- out[, 1:(6 * ntraits), drop = FALSE]
  names(out) <- name_vec
  out <- cbind(crosses_df, out)
  return(out)
}
