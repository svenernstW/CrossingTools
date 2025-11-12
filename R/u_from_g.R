#' Backsolve marker effects (mu) from individual effects (g)
#'
#' Computes marker effects \eqn{\mu = (scalingFactor * M' * G^{-1}) * g}.
#'
#' @param M Numeric matrix (n x p): marker/genotype matrix (individuals in rows, markers in columns).
#' @param G Numeric matrix (n x n): relationship matrix among individuals.
#' @param g Numeric (n x k) or length-n vector: individual effects (one or more traits).
#' @param scaling.factor Numeric scalar used to scale the GRM.
#' @param tol Numeric scalar tolerance for near-singularity checks in \code{G}. Default: \code{1e-10}.
#' @param n.Threads Integer >= 1; OpenMP threads (if enabled at compile time). Default: \code{4L}.
#'
#' @return A numeric \code{p x k} matrix of marker effects (\code{mu_matrix}).
#' @export
u_from_g <- function(M,
                     G,
                     g,
                     scaling.factor,
                     tol = 1e-10,
                     n.Threads = 4L) {
  # Coerce base types
  if (!is.matrix(M)) M <- as.matrix(M)
  if (!is.matrix(G)) G <- as.matrix(G)
  if (is.vector(g))  g <- matrix(as.numeric(g), ncol = 1)
  if (!is.matrix(g)) g <- as.matrix(g)

  # Basic shape checks
  if (!is.numeric(M) || !is.numeric(G) || !is.numeric(g))
    stop("`M`, `G`, and `g` must be numeric.")

  n <- nrow(M); p <- ncol(M)
  if (n < 1L || p < 1L) stop("`M` must have at least 1 row and 1 column.")
  if (nrow(G) != n || ncol(G) != n) stop("`G` must be n x n with n = nrow(M).")
  if (nrow(g) != n) stop("`g` must have nrow(g) == nrow(M).")

  k <- ncol(g)
  if (k < 1L) stop("`g` must have at least one column (trait).")

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
    M                 = M,
    G                 = G,
    g                 = g,
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
