#' Backsolve marker effects from individual effects
#'
#' Backsolves marker effects from individual genetic effects using
#' \deqn{\boldsymbol{\mu} = s\mathbf{M}^{\mathsf{T}}\mathbf{G}^{-1}\mathbf{g},}
#' where \eqn{\mathbf{M}} is the marker matrix, \eqn{\mathbf{G}} is the genomic
#' relationship matrix, \eqn{\mathbf{g}} contains individual effects, and
#' \eqn{s} is a scaling factor determined internally.
#'
#' Multiple traits can be processed simultaneously by supplying
#' \code{genotype.effects} as a matrix or data frame with one column per trait.
#'
#' @param marker.mat Numeric matrix with individuals in rows and markers in
#'   columns.
#' @param G.mat Numeric square matrix containing genomic relationships among
#'   the individuals in \code{marker.mat}.
#' @param genotype.effects Numeric vector, matrix, or data frame containing
#'   individual genetic effects. Rows must correspond to individuals in
#'   \code{marker.mat}; columns represent traits.
#'
#' @return A data frame with markers in rows and traits in columns containing
#'   the backsolved marker effects.
#'
#' @export


backsolve_marker_effects <- function(marker.mat,
                     G.mat,
                     genotype.effects) {
  # Coerce base types
  traits <- names(genotype.effects)
  effects <- as.matrix(genotype.effects)
  G <- G.mat
  tol = 1e-10
  n.Threads = 1

  scaling.factor <- 1

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
  if (nrow(effects) != n) stop("`effects` must have nrow(effects) == nrow(marker.mat).")

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
  res <- cpp_u_from_from_g_simple(
    M                 = marker.mat,
    G                 = G,
    g                 = effects,
    scalingFactor     = scalingFactor
  )

  # Return only the mu matrix
  res <- as.data.frame(res$mu_matrix)

  row.names(res) <- names(marker.mat)

  names(res) <- traits

  res
}
