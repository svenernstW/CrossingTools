#' Calculation of a selection index
#'
#' Calculates a multi-trait selection index from genomic or phenotypic BLUPs.
#' The function supports two alternative index formulations:
#'
#' \itemize{
#'   \item \strong{Desired gains index} , based on user-specified
#'   desired trait gains and a trait (co)variance matrix.
#'   \item \strong{Smith–Hazel index}, based on user-specified economic weights.
#' }
#'
#' Exactly one of \code{gains} or \code{weights} must be provided.
#'
#' @param genotype.effects A numeric matrix (\eqn{n_{genotype} \times n_{trait}}) of multi-trait BLUPs.
#'
#' @param var.mat Optional variance or covariance matrix used for the desired gains index.
#'   If provided, must be one of:
#'   \itemize{
#'     \item A numeric matrix of dimension
#'       \eqn{(n_{genotype} \times n_{trait}) \times (n_{genotype} \times n_{trait})}
#'       containing the posterior variance
#'       \eqn{\mathrm{var}(\tilde{g})} of the stacked BLUPs. In this case, the
#'       \emph{marginal} posterior trait covariance (obtained by averaging
#'       per-genotype trait blocks; Werner et al.) is used for index calculation.
#'     \item A numeric matrix of dimension \eqn{n_{trait} \times n_{trait}}
#'       giving the trait covariance matrix directly.
#'   }
#'   If \code{var.mat} is \code{NULL}, the trait covariance matrix is estimated
#'   empirically as \code{cov(effects)}.
#'
#' @param desired.gains Numeric vector of length \eqn{n_{trait}} specifying the desired gains
#'   for each trait. Used only for the desired gains index.
#'
#' @param weights Numeric vector of length \eqn{n_{trait}} specifying economic
#'   weights for each trait. Used only for the Smith–Hazel index.
#'
#' @param nthreads Integer. Number of OpenMP threads to use.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{index}: A data.frame containing the selection index values for each genotype.
#'   \item \code{weights}: The trait weights used to construct the index.
#' }
#'
#' @export

make_index <- function(
    genotype.effects, weights = NULL, var.mat=NULL, desiredgains = NULL,
    nthreads = 4L
) {
n.Threads <- nthreads
gains <- desired.gains
effects <- genotype.effects
if(!is.null(gains) & !is.null(weights) ){
  stop("Provide either weights or gains, function can only handle one at a time!")
}

if(is.null(gains) & is.null(weights) ){
  stop("Provide either weights or gains!")
}


#  Coerce & basic shapes
if (!is.matrix(effects)) effects <- as.matrix(effects)
nG <- nrow(effects); nT <- ncol(effects)
if (nG < 1L || nT < 1L) stop("`effects` must have at least 1 row and 1 column.")
#Check effects

if (any(!is.finite(effects))) stop("`effects` must contain only finite values.")


# Threads (match parameter name n.Threads)
if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
  stop("`nthreads` must be a positive integer.")
}
nThreads <- as.integer(n.Threads)



if(!is.null(gains) ){

  if (!is.null(var.mat)) {
    var.mat <- as.matrix(var.mat)
    if (!is.numeric(var.mat)) stop("`var.mat` must be numeric.")
    if (nrow(var.mat) != ncol(var.mat)) stop("`var.mat` must be square.")

    if (nrow(var.mat) == nrow(effects)*ncol(effects) && ncol(var.mat) == nrow(effects)*ncol(effects)) {
      use.marginal.V <- TRUE
      use.V.approx   <- FALSE
    }

    if (nrow(var.mat) == ncol(effects) && ncol(var.mat) == ncol(effects)) {
      use.marginal.V <- FALSE
      use.V.approx   <- TRUE
      V.approx <- var.mat
    }

    if(nrow(var.mat)  != nrow(effects)*ncol(effects) & nrow(var.mat) != ncol(effects)){
      stop("`var.mat` must be either nTrait x nTrait or (nG*nT) x (nG*nT).")
    }

  } else {
    use.marginal.V <- FALSE
    use.V.approx   <- TRUE
    V.approx <- cov(effects)
  }










  gains <- as.numeric(gains)
  if (length(gains) != nT || any(!is.finite(gains))) {
    stop("`gains` must be a numeric vector of length ncol(effects) with finite values.")
  }






  #  Call C++
  temp <- cpp_calculate_desired_gains(
    A        = effects,
    V = if (use.marginal.V) var.mat else matrix(0, 1, 1),
    approxV  = if (use.V.approx) V.approx else matrix(0, 1, 1),
    gains    = gains,
    useMargV = use.marginal.V,
    useV     = FALSE,
    useapproxV = use.V.approx,
    nThreads = nThreads
  )
  #  Format return
  out <- list(index = data.frame(index = as.numeric(temp$index)),weights=as.numeric(temp$weight))
  names(out$weights) <- names(effects)
  return(out)

}

if(!is.null(weights) ){
  weights <- as.numeric(weights)
  names(weights) <- names(effects)
  if (length(weights) != nT || any(!is.finite(weights))) {
    stop("`weights` must be a numeric vector of length ncol(effects) with finite values.")
  }
  out <- list(index = data.frame(index = as.numeric(effects %*% weights)),weights=weights)
  return(out)
}


}
