#' Calculation of desired gains index (Werner et al.)
#'
#' Calculates the desired gains index based on multi-trait BLUPs and a chosen
#' (co)variance matrix option.
#'
#' @param effects A matrix (n_genotype x n_trait) of multi-trait BLUPs.
#' @param V Numeric matrix ((n_genotype*n_trait) x (n_genotype*n_trait)) with the
#'   posterior variance var(\eqn{\tilde{g}}) of the provided effects. If provided,
#'   the *marginal* posterior variance (mean of per-genotype trait-blocks from \code{V})
#'   is used to calculate the index. If not provided an empirical estimate of V as cov(A)
#'   is used.
#' @param gains Numeric vector (length = n_trait) of desired gains.
#' @param n.Threads Number of OpenMP threads.
#'
#' @return A list with index values and trait weights.
#' @export
calculate_desired_gains <- function(
    effects, V=NULL, gains = NULL,

    n.Threads = 4L
) {

  if(is.null(V)){
    use.marginal.V = FALSE
    use.V.approx = TRUE
    }
  if(!is.null(V)){
    use.marginal.V = TRUE
    use.V.approx = FALSE
    }

  #  Coerce & basic shapes
  if (!is.matrix(effects)) effects <- as.matrix(effects)
  nG <- nrow(effects); nT <- ncol(effects)
  if (nG < 1L || nT < 1L) stop("`effects` must have at least 1 row and 1 column.")
  #Check effects

  if (any(!is.finite(effects))) stop("`effects` must contain only finite values.")

  # Flags: require exactly one TRUE to avoid ambiguity (C++ overwrites in order)
  flags <- c(use.marginal.V, use.V.approx)
  if (sum(flags) == 0L) stop("Choose exactly one of: use.marginal.V,  use.V.approx.")
  if (sum(flags) > 1L)  stop("Only one of use.marginal.V, use.V.approx can be TRUE.")

  # Gains
  if (is.null(gains)) stop("No gains provided")

  gains <- as.numeric(gains)
  if (length(gains) != nT || any(!is.finite(gains))) {
    stop("`gains` must be a numeric vector of length ncol(effects) with finite values.")
  }

  # V / V.approx checks depending on flag
  if (use.marginal.V) {
    if (missing(V) || is.null(V)) stop("`V` is required when  use.marginal.V is TRUE.")
    if (!is.matrix(V)) V <- as.matrix(V)
    if (!is.numeric(V)) stop("`V` must be numeric.")
    if (nrow(V) != nG * nT || ncol(V) != nG * nT) {
      stop("`V` must be a square matrix of dimension (nrow(effects)*ncol(effects)) x (nrow(effects)*ncol(effects)).")
    }
  }
  if (use.V.approx) {
    V.approx <- cov(effects)

  }

  # Threads (match parameter name n.Threads)
  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

  #  Call C++
  out <- cpp_calculate_desired_gains(
    A        = effects,
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
  rownames(out$index) <- rownames(effects)
  out$weight <- as.numeric(out$weight)
  return(out)
}
