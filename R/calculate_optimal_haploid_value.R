#' Optimal Haplotypic Value (OHV) for proposed crosses
#'
#' Computes the OHV for each cross using block-wise local GEBVs. If \code{haplotype.blocks}
#' is not supplied (or empty), each marker column of \code{marker.mat} is treated as a
#' separate block.
#' @param crosses A numeric matrix or data.frame with two columns (P1, P2), giving
#'   indices of parents (rows of \code{marker.mat}) for each proposed cross.
#' @param haplotype.blocks A list of integer vectors. Each vector contains  marker
#'   indices (columns of \code{marker.mat}) belonging to that block. If \code{NULL} or length 0,
#'   each marker is used as a separate block.
#' @param marker.mat Numeric genotype/marker matrix (individuals x markers).
#' @param effects Numeric matrix of marker effects.
#' @param n.Threads Integer larger 1. OpenMP threads (if enabled at compile time).
#'
#' @return A data.frame with columns \code{OHV}.
#' @export
calculate_optimal_haploid_value <- function(crosses, haplotype.blocks = NULL, marker.mat, effects, n.Threads = 4L) {
  # ---- Normalize inputs ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  if (!is.matrix(crosses)) crosses <- as.matrix(crosses)

  effects <- as.matrix(effects)

  # ---- Basic checks on marker.mat / crosses ----
  if (!is.numeric(marker.mat)) stop("`marker.mat` must be a numeric matrix.")
  if (ncol(marker.mat) < 1L || nrow(marker.mat) < 1L) stop("`marker.mat` must have at least 1 row and 1 column.")

  if (ncol(crosses) != 2L) stop("`crosses` must have exactly 2 columns (P1, P2).")
  if (!all(is.finite(crosses))) stop("`crosses` contains non-finite entries.")
  if (any(crosses < 1 | crosses > nrow(marker.mat))) {
    stop("Some parent indices in `crosses` are outside 1..nrow(marker.mat).")
  }

  # ---- Checks for effects (marker effects matrix) ----
  if (!is.numeric(effects)) stop("`effects` must be numeric (matrix or vector).")
  if (nrow(effects) != ncol(marker.mat)) {
    stop("`effects` must have nrow(effects) == ncol(marker.mat). Got nrow(effects) = ", nrow(effects),
         " but ncol(marker.mat) = ", ncol(marker.mat), ".")
  }
  if (!all(is.finite(effects))) stop("`effects` contains non-finite values.")
  if (ncol(effects) < 1L) stop("`effects` must have at least one column (trait).")

  # ---- haplotype.blocks ----
  if (missing(haplotype.blocks) || is.null(haplotype.blocks) || length(haplotype.blocks) == 0L) {
    haplotype.blocks <- lapply(1:ncol(marker.mat), function(j) j)
  }
  if (!is.list(haplotype.blocks) || length(haplotype.blocks) < 1L) {
    stop("`haplotype.blocks` must be a non-empty list of integer index vectors (or NULL to use single-marker blocks).")
  }
  for (i in seq_along(haplotype.blocks)) {
    blk <- haplotype.blocks[[i]]
    if (length(blk) < 1L) stop("Block ", i, " in `haplotype.blocks` is empty.")
    if (!is.numeric(blk)) stop("Block ", i, " in `haplotype.blocks` must be numeric/integer indices.")
    if (any(!is.finite(blk))) stop("Block ", i, " in `haplotype.blocks` contains non-finite indices.")
    if (any(blk < 1 | blk > ncol(marker.mat))) {
      stop("Block ", i, " in `haplotype.blocks` has indices outside 1..ncol(marker.mat).")
    }
    haplotype.blocks[[i]] <- as.integer(blk)
  }

  # ---- n.Threads ----
  if (length(n.Threads) != 1L || !is.finite(n.Threads) ||
      abs(n.Threads - round(n.Threads)) > .Machine$double.eps^0.5 || n.Threads < 1) {
    stop("`n.threads` must be a single integer >= 1.")
  }
  n.Threads <- as.integer(n.Threads)

  # ---- Map to C++ parameter names (do NOT change C++) ----
  Crosses  <- crosses
  HBlocks  <- haplotype.blocks
  nThreads <- n.Threads
  temp <- list()

  for (i in 1:ncol(effects)) {
    temp[[i]] <- cpp_calcOHV(
      Crosses  = Crosses,
      HBlocks  = HBlocks,
      M        = marker.mat,
      mu_vec   = effects[, i],
      nThreads = nThreads
    )
  }

  out <- as.data.frame(do.call(cbind, temp))

  names(out) <- paste0("OHV_", 1:ncol(effects))
  out <- cbind(crosses, out)

  return(out)
}
