#' Calculation of desired gains index (Werner et al.)
#'
#' Calculates the desired gains index based on multi-trait BLUPs and a chosen
#' (co)variance matrix option.
#'
#' @param A A matrix (n_genotype x n_trait) of multi-trait BLUPs.
#' @param V Numeric matrix ((n_genotype*n_trait) x (n_genotype*n_trait)) with the
#'   posterior variance var(\eqn{\tilde{g}}).
#'   Required if \code{use.marginal.V} is TRUE.
#' @param gains Numeric vector (length = n_trait) of desired gains.
#' @param use.marginal.V Logical. If TRUE, uses the *marginal* posterior variance
#'   (mean of per-genotype trait-blocks from \code{V}).
#' @param use.V.approx Logical. If TRUE, uses an empirical estimate of V as cov(A).
#' @param n.Threads Number of OpenMP threads.
#'
#' @return A list with index values and trait weights.
#' @export
calculate_desired_gains <- function(
    A, V, gains = NULL,
    use.marginal.V = FALSE, use.V.approx = FALSE,
    n.Threads = 4L
) {
  #  Coerce & basic shapes
  if (!is.matrix(A)) A <- as.matrix(A)
  nG <- nrow(A); nT <- ncol(A)
  if (nG < 1L || nT < 1L) stop("`A` must have at least 1 row and 1 column.")
  #Check A

  if (any(!is.finite(A))) stop("`A` must contain only finite values.")

  # Flags: require exactly one TRUE to avoid ambiguity (C++ overwrites in order)
  flags <- c(use.marginal.V, use.V.approx)
  if (sum(flags) == 0L) stop("Choose exactly one of: use.marginal.V,  use.V.approx.")
  if (sum(flags) > 1L)  stop("Only one of use.marginal.V, use.V.approx can be TRUE.")

  # Gains
  if (is.null(gains)) stop("No gains provided")

  gains <- as.numeric(gains)
  if (length(gains) != nT || any(!is.finite(gains))) {
    stop("`gains` must be a numeric vector of length ncol(A) with finite values.")
  }

  # V / V.approx checks depending on flag
  if (use.marginal.V) {
    if (missing(V) || is.null(V)) stop("`V` is required when  use.marginal.V is TRUE.")
    if (!is.matrix(V)) V <- as.matrix(V)
    if (!is.numeric(V)) stop("`V` must be numeric.")
    if (nrow(V) != nG * nT || ncol(V) != nG * nT) {
      stop("`V` must be a square matrix of dimension (nrow(A)*ncol(A)) x (nrow(A)*ncol(A)).")
    }
  }
  if (use.V.approx) {
    V.approx <- cov(A)

  }

  # Threads (match parameter name n.Threads)
  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  #  Call C++
  out <- cpp_calculate_desired_gains(
    A        = A,
    V = if (use.marginal.V) V else matrix(0, 1, 1),
    approxV  = if (use.V.approx) V.approx else matrix(0, 1, 1),
    gains    = gains,
    useMargV = use.marginal.V,
    useV     = FALSE,
    useapproxV = use.V.approx,
    nThreads = nThreads
  )

  #  Format return
  out$index <- data.frame(DG_index = as.numeric(out$index))
  rownames(out$index) <- rownames(A)
  out$weight <- as.numeric(out$weight)
  return(out)
}
