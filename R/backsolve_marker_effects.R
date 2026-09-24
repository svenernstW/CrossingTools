#' Backsolve marker effects from individual genetic effects
#'
#' Backsolves marker effects from predicted individual genetic effects using
#' the equivalence between genomic relationship-matrix and marker-effect models
#' (Stranden and Garrick, 2009).
#'
#' Let \eqn{M} denote the marker design matrix, \eqn{G} the genomic relationship
#' matrix, and \eqn{g_t} the vector of predicted individual genetic effects for
#' trait \eqn{t}. Marker effects are first obtained as
#'
#' \deqn{
#' \tilde{\mu}_t = M^\top G^{-1} g_t .
#' }
#'
#' Because genomic relationship matrices may differ from \eqn{M M^\top} by a
#' constant scaling factor, the backsolved marker effects are subsequently
#' rescaled so that the genetic values reconstructed from the marker effects
#' optimally reproduce the supplied individual genetic effects. For each trait,
#' the scalar
#'
#' \deqn{
#' c_t =
#' \frac{(M\tilde{\mu}_t)^\top g_t}
#'      {(M\tilde{\mu}_t)^\top(M\tilde{\mu}_t)}
#' }
#'
#' is calculated, and the final marker effects are
#'
#' \deqn{
#' \mu_t = c_t \tilde{\mu}_t .
#' }
#'
#' Thus, scaling is estimated independently for each trait and corresponds to a
#' single constant applied to all marker effects within that trait. This is
#' appropriate when \code{G.mat} is proportional to the crossproduct of the
#' supplied marker matrix, for example
#'
#' \deqn{
#' G = \frac{M M^\top}{s},
#' }
#'
#' where \eqn{s} is a scalar normalization constant.
#'
#' Genomic relationship matrices using marker-specific or other non-scalar
#' scaling are not supported.
#'
#' The marker matrix supplied to the function must use the same coding,
#' centring, and marker representation as that used to construct
#' \code{G.mat}. For example, centred additive marker genotypes can be used to
#' recover average allele-substitution effects from genomic breeding values
#' when \code{G.mat} was constructed from the same centred additive marker
#' matrix (VanRaden, 2008; Stranden and Garrick, 2009). Likewise, dominance
#' effects can be backsolved by supplying the corresponding dominance design
#' matrix and dominance genomic relationship matrix.
#'
#' Multiple traits can be processed simultaneously by supplying
#' \code{genotype.effects} as a matrix or data frame with one column per trait.
#' Marker effects and their scalar rescaling factors are calculated separately
#' for each trait using the same marker and genomic relationship matrices.
#'
#' @param marker.mat Numeric marker design matrix with individuals in rows and
#'   markers in columns. Its coding, centring, and marker representation must
#'   correspond to those used to construct \code{G.mat}.
#' @param G.mat Numeric square genomic relationship matrix among the individuals
#'   in \code{marker.mat}.
#' @param genotype.effects Numeric vector, matrix, or data frame containing
#'   predicted individual genetic effects. Rows must correspond to individuals
#'   in \code{marker.mat}; columns represent traits.
#'
#' @return A data frame with markers in rows and traits in columns containing
#'   the backsolved marker effects.
#'
#' @references
#' VanRaden, P. M. (2008).
#' Efficient methods to compute genomic predictions.
#' \emph{Journal of Dairy Science}, 91(11), 4414--4423.
#' \doi{10.3168/jds.2007-0980}
#'
#' Stranden, I. and Garrick, D. J. (2009).
#' Derivation of equivalent computing algorithms for genomic predictions and
#' reliabilities of animal merit.
#' \emph{Journal of Dairy Science}, 92(6), 2971--2975.
#' \doi{10.3168/jds.2008-1929}
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
