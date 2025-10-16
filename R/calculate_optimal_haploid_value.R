#' Optimal Haplotypic Value (OHV) for proposed crosses
#'
#' Computes the OHV for each cross using block-wise local GEBVs. If \code{haplotype.blocks}
#' is not supplied (or empty), each marker column of \code{M} is treated as a
#' separate block.
#' @param crosses A numeric matrix or data.frame with two columns (P1, P2), giving
#'   indices of parents (rows of \code{M}) for each proposed cross.
#' @param haplotype.blocks A list of integer vectors. Each vector contains  marker
#'   indices (columns of \code{M}) belonging to that block. If \code{NULL} or length 0,
#'   each marker is used as a separate block.
#' @param M Numeric genotype/marker matrix (individuals x markers).
#' @param U Numeric vector of marker effects (length \code{ncol(M)}).
#' @param nthreads Integer \(\ge 1\). OpenMP threads (if enabled at compile time).
#'
#' @return A data.frame with column \code{OHV}.
#' @export
calculate_optimal_haploid_value <- function(crosses, haplotype.blocks=NULL, M, U, nthreads = 4L) {
  # ---- Normalize inputs ----
  if (!is.matrix(M)) M <- as.matrix(M)
  if (!is.matrix(crosses)) crosses <- as.matrix(crosses)

  U <- as.matrix(U)

  # ---- Basic checks on M / crosses ----
  if (!is.numeric(M)) stop("`M` must be a numeric matrix.")
  if (ncol(M) < 1L || nrow(M) < 1L) stop("`M` must have at least 1 row and 1 column.")

  if (ncol(crosses) != 2L) stop("`crosses` must have exactly 2 columns (P1, P2).")
  if (!all(is.finite(crosses))) stop("`crosses` contains non-finite entries.")
  if (any(crosses < 1 | crosses > nrow(M))) {
    stop("Some parent indices in `crosses` are outside 1..nrow(M).")
  }

  # ---- Checks for U (marker effects matrix) ----
  # Expect one row per marker in M and one column per trait/index
  if (!is.numeric(U)) stop("`U` must be numeric (matrix or vector).")
  if (nrow(U) != ncol(M)) {
    stop("`U` must have nrow(U) == ncol(M). Got nrow(U) = ", nrow(U),
         " but ncol(M) = ", ncol(M), ".")
  }
  if (!all(is.finite(U))) stop("`U` contains non-finite values.")
  if (ncol(U) < 1L) stop("`U` must have at least one column (trait).")

  # ---- haplotype.blocks ----
  if (missing(haplotype.blocks) || is.null(haplotype.blocks) || length(haplotype.blocks) == 0L) {
    haplotype.blocks <- lapply(1:ncol(M), function(j) j)
  }
  if (!is.list(haplotype.blocks) || length(haplotype.blocks) < 1L) {
    stop("`haplotype.blocks` must be a non-empty list of integer index vectors (or NULL to use single-marker blocks).")
  }
  for (i in seq_along(haplotype.blocks)) {
    blk <- haplotype.blocks[[i]]
    if (length(blk) < 1L) stop("Block ", i, " in `haplotype.blocks` is empty.")
    if (!is.numeric(blk)) stop("Block ", i, " in `haplotype.blocks` must be numeric/integer indices.")
    if (any(!is.finite(blk))) stop("Block ", i, " in `haplotype.blocks` contains non-finite indices.")
    if (any(blk < 1 | blk > ncol(M))) {
      stop("Block ", i, " in `haplotype.blocks` has indices outside 1..ncol(M).")
    }
    # make sure indices are integers
    haplotype.blocks[[i]] <- as.integer(blk)
  }

  # ---- nthreads ----
  if (length(nthreads) != 1L || !is.finite(nthreads) ||
      abs(nthreads - round(nthreads)) > .Machine$double.eps^0.5 || nthreads < 1) {
    stop("`nthreads` must be a single integer >= 1.")
  }
  nthreads <- as.integer(nthreads)

  # ---- Map to C++ parameter names (do NOT change C++) ----
  Crosses  <- crosses
  HBlocks  <- haplotype.blocks
  nThreads <- nthreads
  temp <- list()
  # ---- Call C++ ----
  for(i in 1:ncol(U)){
    temp[[i]] <- cpp_calcOHV(
      Crosses  = Crosses,
      HBlocks  = HBlocks,
      M        = M,
      mu_vec        = U[,i],
      nThreads = nThreads
    )
  }



  out <- as.data.frame(do.call(cbind,temp))

  names(out) <- paste0("OHV_", 1:ncol(U))


  return(out)
}
