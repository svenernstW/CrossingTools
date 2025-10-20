#' Backsolve marker effects (mu) from individual effects (g)
#'
#' Computes marker effects \eqn{\mu = (scalingFactor * M' * G^{-1}) * g}, with
#' optional prior/posterior variance-covariance outputs and LD-variance by groups.
#'
#' @param M Numeric matrix (n x p): marker/genotype matrix (individuals in rows, markers in columns).
#' @param G Numeric matrix (n x n): relationship matrix among individuals.
#' @param g Numeric (n x k) or length-n vector: individual effects (one or more traits).
#' @param scalingFactor Numeric scalar multiplier.
#' @param tol Numeric scalar tolerance for near-singularity checks in \code{G}.
#' @param LDvar Logical; if \code{TRUE}, compute LD variance by groups.
#' @param Grouping \code{NULL} or a data.frame with a column \code{chr} of length \code{ncol(M)}
#'   giving group membership (only used if \code{LDvar=TRUE}). If \code{NULL} with \code{LDvar=TRUE},
#'   all markers are treated as one group. If supplied be aware that the order must match exaclty that of M
#' @param calcPriorVcov Logical; if \code{TRUE}, compute prior Vcov over marker-by-trait effects.
#' @param PriorVcov \code{NULL}, \code{p x p} or \code{(p*k) x (p*k)} numeric matrix; see C++ docstring.
#' @param sigmasq \code{NULL}, scalar, length-k vector, or \code{k x k} matrix specifying trait co-variance(s);
#'   required if \code{calcPriorVcov=TRUE}.
#' @param calcPosteriorVcov Logical; if \code{TRUE}, compute posterior Vcov using \code{PEV}.
#' @param PEV \code{NULL} or numeric matrix of size \code{(n*k) x (n*k)} (prediction error variance).
#' @param nThreads Integer larger 1; OpenMP threads (if enabled at compile time).
#'
#' @return A list with elements matching the C++ return:
#' \describe{
#'   \item{\code{mu_matrix}}{\code{p x k} matrix of marker effects.}
#'   \item{\code{variance_LD}}{LD variance by groups (matrix if \code{k=1}, else list), or \code{NULL} if \code{LDvar=FALSE}.}
#'   \item{\code{prior_Vcov}}{Prior Vcov matrix if requested, else \code{NULL}.}
#'   \item{\code{posterior_Vcov}}{Posterior Vcov matrix if requested, else empty \code{0x0} matrix.}
#'   \item{\code{group_levels}}{Group labels used (only when \code{LDvar=TRUE}).}
#' }
#' @export
u_from_g <- function(M,
                          G,
                          g,
                          scalingFactor,
                          tol = 1e-10,
                          LDvar = FALSE,
                          Grouping = NULL,
                          calcPriorVcov = FALSE,
                          PriorVcov = NULL,
                          sigmasq = NULL,
                          calcPosteriorVcov = FALSE,
                          PEV = NULL,
                          nThreads = 4L) {
  #  Coerce base types (keep your arg names)
  if (!is.matrix(M)) M <- as.matrix(M)
  if (!is.matrix(G)) G <- as.matrix(G)
  if (is.vector(g))  g <- matrix(as.numeric(g), ncol = 1)
  if (!is.matrix(g)) g <- as.matrix(g)

  #  Basic shape checks
  if (!is.numeric(M) || !is.numeric(G) || !is.numeric(g))
    stop("`M`, `G`, and `g` must be numeric.")

  n <- nrow(M); p <- ncol(M)
  if (n < 1L || p < 1L) stop("`M` must have at least 1 row and 1 column.")
  if (nrow(G) != n || ncol(G) != n) stop("`G` must be n x n with n = nrow(M).")
  if (nrow(g) != n) stop("`g` must have nrow(g) == nrow(M).")

  k <- ncol(g)
  if (k < 1L) stop("`g` must have at least one column (trait).")

  # Scalars
  if (length(scalingFactor) != 1L || !is.finite(scalingFactor))
    stop("`scalingFactor` must be a single finite number.")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("`tol` must be a single positive number.")
  if (length(LDvar) != 1L || !is.logical(LDvar))
    stop("`LDvar` must be a single logical.")
  if (length(calcPriorVcov) != 1L || !is.logical(calcPriorVcov))
    stop("`calcPriorVcov` must be a single logical.")
  if (length(calcPosteriorVcov) != 1L || !is.logical(calcPosteriorVcov))
    stop("`calcPosteriorVcov` must be a single logical.")
  if (length(nThreads) != 1L || !is.finite(nThreads) ||
      abs(nThreads - round(nThreads)) > .Machine$double.eps^0.5 || nThreads < 1L)
    stop("`nThreads` must be a single integer >= 1.")
  nThreads <- as.integer(nThreads)

  #  Grouping (used only when LDvar=TRUE)
  if (isTRUE(LDvar)) {
    if (is.null(Grouping)) {
      # default: single group
      Grouping <- data.frame(chr = rep(1L, p))
    } else if (is.vector(Grouping) && !is.list(Grouping)) {
      if (length(Grouping) != p) stop("`Grouping` vector must have length ncol(M).")
      Grouping <- data.frame(chr = as.integer(Grouping))
    } else if (is.data.frame(Grouping)) {
      if (!("chr" %in% names(Grouping)))
        stop("`Grouping` data.frame must contain a 'chr' column.")
      if (nrow(Grouping) != p)
        stop("`Grouping` must have one row per marker (nrow(Grouping) == ncol(M)).")
      Grouping$chr <- as.integer(Grouping$chr)
    } else {
      stop("`Grouping` must be NULL, a vector of group labels, or a data.frame with column 'chr'.")
    }
  } else {
    # When LDvar=FALSE, we pass NULL through.
    Grouping <- NULL
  }

  #  Prior Vcov requirements
  if (isTRUE(calcPriorVcov) && is.null(sigmasq)) {
    stop("`calcPriorVcov = TRUE` requires `sigmasq` (scalar, length-k vector, or k x k matrix).")
  }
  # If PriorVcov provided, sanity-check its dims (p x p) or (p*k x p*k)
  if (!is.null(PriorVcov)) {
    PriorVcov <- as.matrix(PriorVcov)
    if (!((nrow(PriorVcov) == p && ncol(PriorVcov) == p) ||
          (nrow(PriorVcov) == p*k && ncol(PriorVcov) == p*k))) {
      stop("`PriorVcov` must be either p x p or (p*k) x (p*k).")
    }
  }

  #  Posterior Vcov requirements
  if (isTRUE(calcPosteriorVcov)) {
    if (is.null(PEV)) stop("`calcPosteriorVcov = TRUE` requires `PEV`.")
    PEV <- as.matrix(PEV)
    if (!(nrow(PEV) == n*k && ncol(PEV) == n*k)) {
      stop("`PEV` must be of size (n*k) x (n*k).")
    }
  }

  #  Call C++
  res <- cpp_u_from_from_g(
    M                  = M,
    G                  = G,
    g                  = g,
    scalingFactor      = scalingFactor,
    tol                = tol,
    LDvar              = LDvar,
    Grouping           = if (is.null(Grouping)) NULL else Grouping,
    calcPriorVcov      = calcPriorVcov,
    PriorVcov          = if (is.null(PriorVcov)) NULL else PriorVcov,
    sigmasq            = if (is.null(sigmasq)) NULL else sigmasq,
    calcPosteriorVcov  = calcPosteriorVcov,
    PEV                = if (is.null(PEV)) NULL else PEV,
    nThreads           = nThreads
  )

  res
}
