#' Backsolve marker effects (mu) from individual effects (g)
#'
#' Computes marker effects \eqn{\mu = (scalingFactor * M' * G^{-1}) * g}, with
#' optional prior/posterior variance-covariance outputs and LD-variance by groups.
#'
#' @param M Numeric matrix (n x p): marker/genotype matrix (individuals in rows, markers in columns).
#' @param G Numeric matrix (n x n): relationship matrix among individuals.
#' @param g Numeric (n x k) or length-n vector: individual effects (one or more traits).
#' @param scaling.factor Numeric scalar, that is used to scale the GRM.
#' @param tol Numeric scalar tolerance for near-singularity checks in \code{G}.
#' @param LD.var Logical; if \code{TRUE}, compute LD variance by groups.
#' @param grouping \code{NULL} or a data.frame with a column \code{chr} of length \code{ncol(M)}
#'   giving group membership (only used if \code{LD.var=TRUE}). If \code{NULL} with \code{LD.var=TRUE},
#'   all markers are treated as one group. If supplied be aware that the order must match exaclty that of M
#' @param calc.prior.Vcov Logical; if \code{TRUE}, compute prior Vcov over marker-by-trait effects.
#' @param prior.Vcov \code{NULL}, \code{p x p} or \code{(p*k) x (p*k)} numeric matrix.
#' @param sigma.sq \code{NULL}, scalar, length-k vector, or \code{k x k} matrix specifying trait co-variance(s);
#'   required if \code{calc.prior.Vcov=TRUE}.
#' @param calc.posterior.Vcov Logical; if \code{TRUE}, compute posterior Vcov using \code{PEV}.
#' @param PEV \code{NULL} or numeric matrix of size \code{(n*k) x (n*k)} (prediction error variance).
#' @param n.Threads Integer larger 1; OpenMP threads (if enabled at compile time).
#'
#' @return A list with elements matching the C++ return:
#' \describe{
#'   \item{\code{mu_matrix}}{\code{p x k} matrix of marker effects.}
#'   \item{\code{variance_LD}}{LD variance by groups (matrix if \code{k=1}, else list), or \code{NULL} if \code{LD.var=FALSE}.}
#'   \item{\code{prior_Vcov}}{Prior Vcov matrix if requested, else \code{NULL}.}
#'   \item{\code{posterior_Vcov}}{Posterior Vcov matrix if requested, else empty \code{0x0} matrix.}
#'   \item{\code{group_levels}}{Group labels used (only when \code{LD.var=TRUE}).}
#' }
#' @export
u_from_g <- function(M,
                     G,
                     g,
                     scaling.factor,
                     tol = 1e-10,
                     LD.var = FALSE,
                     grouping = NULL,
                     calc.prior.Vcov = FALSE,
                     prior.Vcov = NULL,
                     sigma.sq = NULL,
                     calc.posterior.Vcov = FALSE,
                     PEV = NULL,
                     n.Threads = 4L) {
  #  Coerce base types
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

  # Scalars (match new arg names; create internal variables expected by C++)
  if (length(scaling.factor) != 1L || !is.finite(scaling.factor))
    stop("`scaling.factor` must be a single finite number.")
  scalingFactor <- as.numeric(scaling.factor)

  if (length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("`tol` must be a single positive number.")

  if (length(LD.var) != 1L || !is.logical(LD.var))
    stop("`LD.var` must be a single logical.")
  LDvar <- LD.var

  if (length(calc.prior.Vcov) != 1L || !is.logical(calc.prior.Vcov))
    stop("`calc.prior.Vcov` must be a single logical.")
  calcPriorVcov <- calc.prior.Vcov

  if (length(calc.posterior.Vcov) != 1L || !is.logical(calc.posterior.Vcov))
    stop("`calc.posterior.Vcov` must be a single logical.")
  calcPosteriorVcov <- calc.posterior.Vcov

  if (length(n.Threads) != 1L || !is.finite(n.Threads) ||
      abs(n.Threads - round(n.Threads)) > .Machine$double.eps^0.5 || n.Threads < 1L)
    stop("`n.Threads` must be a single integer >= 1.")
  nThreads <- as.integer(n.Threads)

  #  Grouping (used only when LDvar=TRUE)
  if (isTRUE(LDvar)) {
    if (is.null(grouping)) {
      # default: single group
      Grouping <- data.frame(chr = rep(1L, p))
    } else if (is.vector(grouping) && !is.list(grouping)) {
      if (length(grouping) != p) stop("`grouping` vector must have length ncol(M).")
      Grouping <- data.frame(chr = as.integer(grouping))
    } else if (is.data.frame(grouping)) {
      if (!("chr" %in% names(grouping)))
        stop("`grouping` data.frame must contain a 'chr' column.")
      if (nrow(grouping) != p)
        stop("`grouping` must have one row per marker (nrow(grouping) == ncol(M)).")
      Grouping <- grouping
      Grouping$chr <- as.integer(Grouping$chr)
    } else {
      stop("`grouping` must be NULL, a vector of group labels, or a data.frame with column 'chr'.")
    }
  } else {
    Grouping <- NULL
  }

  #  Prior Vcov requirements
  sigmasq <- sigma.sq
  if (isTRUE(calcPriorVcov) && is.null(sigmasq)) {
    stop("`calc.prior.Vcov = TRUE` requires `sigma.sq` (scalar, length-k vector, or k x k matrix).")
  }
  PriorVcov <- prior.Vcov
  if (!is.null(PriorVcov)) {
    PriorVcov <- as.matrix(PriorVcov)
    if (!((nrow(PriorVcov) == p && ncol(PriorVcov) == p) ||
          (nrow(PriorVcov) == p*k && ncol(PriorVcov) == p*k))) {
      stop("`prior.Vcov` must be either p x p or (p*k) x (p*k).")
    }
  }

  #  Posterior Vcov requirements
  if (isTRUE(calcPosteriorVcov)) {
    if (is.null(PEV)) stop("`calc.posterior.Vcov = TRUE` requires `PEV`.")
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
