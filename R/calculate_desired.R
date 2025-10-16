#' Calculation of desired gains index (Werner et al.)
#'
#' Calculates the desired gains index based on multi-trait BLUPs and a chosen
#' (co)variance matrix option.
#'
#' @param A A matrix (n_genotype x n_trait) of multi-trait BLUPs.
#' @param V Numeric matrix ((n_genotype*n_trait) x (n_genotype*n_trait)) with the
#'   posterior variance var(\eqn{\tilde{g}}).
#'   Required if \code{use.V} or \code{use.marginal.V} is TRUE.
#' @param V.approx Numeric matrix (n_trait x n_trait) giving an approximation of
#'   var(\eqn{\tilde{g}}), e.g. \code{cov(A)}. Required if \code{use.V.approx} is TRUE.
#' @param gains Numeric vector (length = n_trait) of desired gains. If NULL, uses
#'   equal weights (\code{rep(1, n_trait)}).
#' @param use.marginal.V Logical. If TRUE, uses the *marginal* posterior variance
#'   (mean of per-genotype trait-blocks from \code{V}).
#' @param use.V Logical. If TRUE, uses each genotype's own trait-block \code{V_g}
#'   to compute a per-genotype index (i.e., \eqn{w_g = V_g^{-1} d}).
#' @param use.V.approx Logical. If TRUE, uses \code{V.approx}.
#' @param nThreads Number of OpenMP threads.
#'
#' @return A one-column data.frame with the desired-gains index for each genotype.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' nG <- 50; nT <- 3
#' A <- matrix(rnorm(nG * nT), nG, nT)               # BLUPs
#' V.approx <- cov(A)                                 # simple approximation
#' V <- kronecker(diag(nG), V.approx)                 # block-diagonal (nG*nT x nG*nT)
#' gains <- c(1, 0.5, 0.2)
#'
#' # Using the approximation:
#' idx1 <- calculate_desired_gains(A, V, V.approx, gains,
#'                                 use.V = FALSE, use.marginal.V = FALSE,
#'                                 use.V.approx = TRUE)
#'
#' # Using the marginal V (averaged across genotype blocks):
#' idx2 <- calculate_desired_gains(A, V, V.approx, gains,
#'                                 use.V = FALSE, use.marginal.V = TRUE,
#'                                 use.V.approx = FALSE)
#'
#' # Using per-genotype V blocks:
#' idx3 <- calculate_desired_gains(A, V, V.approx, gains,
#'                                 use.V = TRUE, use.marginal.V = FALSE,
#'                                 use.V.approx = FALSE)
#' }
#' @export
calculate_desired_gains <- function(
    A, V, V.approx, gains = NULL,
    use.V = TRUE, use.marginal.V = FALSE, use.V.approx = FALSE,
    nThreads = 4L
) {
  #  Coerce & basic shapes 
  if (!is.matrix(A)) A <- as.matrix(A)
  nG <- nrow(A); nT <- ncol(A)
  if (nG < 1L || nT < 1L) stop("`A` must have at least 1 row and 1 column.")
  
  # Flags: require exactly one TRUE to avoid ambiguity (C++ overwrites in order)
  flags <- c(use.marginal.V, use.V, use.V.approx)
  if (sum(flags) == 0L) stop("Choose exactly one of: use.marginal.V, use.V, use.V.approx.")
  if (sum(flags) > 1L)  stop("Only one of use.marginal.V, use.V, use.V.approx can be TRUE.")
  
  # Gains
  if (is.null(gains)) gains <- rep(1, nT)
  gains <- as.numeric(gains)
  if (length(gains) != nT || any(!is.finite(gains))) {
    stop("`gains` must be a numeric vector of length ncol(A) with finite values.")
  }
  
  # V / V.approx checks depending on flag
  if (use.V || use.marginal.V) {
    if (missing(V) || is.null(V)) stop("`V` is required when use.V or use.marginal.V is TRUE.")
    if (!is.matrix(V)) V <- as.matrix(V)
    if (!is.numeric(V)) stop("`V` must be numeric.")
    if (nrow(V) != nG * nT || ncol(V) != nG * nT) {
      stop("`V` must be a square matrix of dimension (nrow(A)*ncol(A)) x (nrow(A)*ncol(A)).")
    }
  }
  if (use.V.approx) {
    if (missing(V.approx) || is.null(V.approx)) stop("`V.approx` is required when use.V.approx is TRUE.")
    if (!is.matrix(V.approx)) V.approx <- as.matrix(V.approx)
    if (!is.numeric(V.approx)) stop("`V.approx` must be numeric.")
    if (!all(dim(V.approx) == c(nT, nT))) {
      stop("`V.approx` must be an ncol(A) x ncol(A) matrix.")
    }
  }
  
  # Threads
  if (length(nThreads) != 1L || !is.finite(nThreads) || nThreads < 1 || nThreads != as.integer(nThreads)) {
    stop("`nThreads` must be a positive integer.")
  }
  nThreads <- as.integer(nThreads)
  
  #  Call C++ 
  temp <- cpp_calculate_desired_gains(
    A        = A,
    V        = if (use.V || use.marginal.V) V else matrix(0, 1, 1),
    approxV  = if (use.V.approx) V.approx else matrix(0, 1, 1),
    gains    = gains,
    useMargV = use.marginal.V,
    useV     = use.V,
    useapproxV = use.V.approx,
    nThreads = nThreads
  )
  
  #  Format return 
  out <- data.frame(DG_index = as.numeric(temp))
  rownames(out) <- rownames(A)
  return(out)
}
