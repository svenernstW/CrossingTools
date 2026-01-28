#' Backsolve marker effects (mu) from individual effects (g)
#'
#' Computes marker effects \eqn{\mu = (scalingFactor * M' * G^{-1}) * g}.
#'
#' @param marker.mat Numeric matrix (n x p): marker/genotype matrix (individuals in rows, markers in columns).
#' @param G.mat Numeric matrix (n x n): relationship matrix among individuals.
#' @param effects Numeric (n x k) or length-n vector: individual effects (one or more traits).
#' @param scaling.factor Numeric scalar used to scale the GRM.
#' @param n.Threads Integer >= 1; OpenMP threads (if enabled at compile time). Default: \code{4L}.
#' @return A numeric \code{p x k} matrix of marker effects (\code{mu_matrix}).
#' @export
get_marker_effects <- function(marker.mat,
                     G.mat,
                     effects,
                     scaling.factor,

                     n.Threads = 4L) {
  # Coerce base types
  G <- G.mat
  tol = 1e-10
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  if (!is.matrix(G)) G <- as.matrix(G)
  if (is.vector(effects))  effects <- matrix(as.numeric(effects), ncol = 1)
  if (!is.matrix(effects)) effects <- as.matrix(effects)

  # Basic shape checks
  if (!is.numeric(marker.mat) || !is.numeric(G) || !is.numeric(effects))
    stop("`marker.mat`, `G`, and `effects` must be numeric.")

  n <- nrow(marker.mat); p <- ncol(marker.mat)
  if (n < 1L || p < 1L) stop("`marker.mat` must have at least 1 row and 1 column.")
  if (nrow(G) != n || ncol(G) != n) stop("`G` must be n x n with n = nrow(marker.mat).")
  if (nrow(marker.mat) != n) stop("`effects` must have nrow(effects) == nrow(marker.mat).")

  k <- ncol(effects)
  if (k < 1L) stop("`effects` must have at least one column (trait).")

  # Scalars
  if (length(scaling.factor) != 1L || !is.finite(scaling.factor))
    stop("`scaling.factor` must be a single finite number.")
  scalingFactor <- as.numeric(scaling.factor)

  if (length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("`tol` must be a single positive number.")

  if (length(n.Threads) != 1L || !is.finite(n.Threads) ||
      abs(n.Threads - round(n.Threads)) > .Machine$double.eps^0.5 || n.Threads < 1L)
    stop("`n.Threads` must be a single integer >= 1.")
  nThreads <- as.integer(n.Threads)

  # Force all optional computations OFF inside the C++ call
  res <- cpp_u_from_from_g(
    M                 = marker.mat,
    G                 = G,
    g                 = effects,
    scalingFactor     = scalingFactor,
    tol               = tol,
    LDvar             = FALSE,
    Grouping          = NULL,
    calcPriorVcov     = FALSE,
    PriorVcov         = NULL,
    sigmasq           = NULL,
    calcPosteriorVcov = FALSE,
    PEV               = NULL,
    nThreads          = nThreads
  )

  # Return only the mu matrix
  res$mu_matrix
}
