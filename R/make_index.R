#' Calculate a multi-trait selection index
#'
#' Calculates a selection index from multi-trait genomic predicted effects
#' using either a desired-gains index or a Smith-Hazel index.
#'
#' For the desired-gains index, trait weights are derived from the specified
#' desired gains and a trait covariance matrix. For the Smith-Hazel index, the
#' supplied trait weights are used directly. Exactly one of
#' \code{desired.gains} or \code{weights} must be provided.
#'
#' @param genotype.effects Numeric matrix or data frame containing genotype
#'   effects, with genotypes in rows and traits in columns.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, a Smith-Hazel index is calculated as a weighted sum of the trait
#'   effects.
#' @param var.mat Optional covariance matrix used for the desired-gains index.
#'   It may be either a trait covariance matrix or the posterior covariance
#'   matrix of the stacked genotype effects. If \code{NULL}, the trait
#'   covariance matrix is estimated from \code{genotype.effects}.
#' @param desired.gains Optional numeric vector specifying the desired gain for
#'   each trait.
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{index}}{A data frame containing the selection index value for
#'     each genotype.}
#'     \item{\code{weights}}{The trait weights used to construct the index.}
#'   }
#'
#' @export


make_index <- function(
    genotype.effects, weights = NULL, var.mat=NULL, desired.gains = NULL,
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
