#' Calculate a multi-trait selection index
#'
#' Constructs a multi-trait selection index from genomic predicted effects
#' using either a Desired Gains Index or user-supplied coefficients for a
#' Smith--Hazel index.
#'
#' Linear selection indices combine information from multiple traits into a
#' single criterion on which individuals can be ranked and selected
#' (Smith, 1936; Hazel, 1943). When \code{weights} are supplied, the index for
#' individual \eqn{i} is calculated as
#' \deqn{I_i = b^\top g_i,}
#' where \eqn{g_i} is the vector of predicted effects for individual \eqn{i}
#' and \eqn{b} contains the supplied index coefficients. This corresponds to
#' the linear-index formulation underlying the Smith--Hazel selection index.
#' The coefficients are supplied directly and are not estimated from economic
#' weights and separate phenotypic and genetic covariance matrices.
#'
#' Alternatively, a Desired Gains Index can be constructed by specifying the
#' desired relative response for each trait. The Desired Gains Index provides
#' a way to express the breeding objective directly in terms of the relative
#' changes desired among traits, rather than specifying index coefficients
#' directly (Werner et al., 2024).
#'
#' When \code{desired.gains} is supplied, the index coefficients are derived
#' from the desired-gain vector and the covariance structure of the predicted
#' genotype effects. The covariance information can be supplied either as a
#' trait covariance matrix or as the covariance matrix of the stacked
#' multi-trait genotype effects. If \code{var.mat = NULL}, the trait covariance
#' matrix is estimated from \code{genotype.effects}.
#'
#' Exactly one of \code{desired.gains} or \code{weights} must be supplied.
#'
#' @param genotype.effects Numeric matrix or data frame containing predicted
#'   genotype effects, with genotypes in rows and traits in columns.
#' @param weights Optional numeric vector with one index coefficient per trait.
#'   When supplied, a Smith--Hazel-type linear selection index is calculated as
#'   the weighted sum of the trait effects. The supplied values are interpreted
#'   directly as index coefficients.
#' @param var.mat Optional numeric covariance matrix used to construct the
#'   Desired Gains Index. It may be either an
#'   \eqn{n_{\mathrm{trait}} \times n_{\mathrm{trait}}} trait covariance
#'   matrix or the covariance matrix of the stacked multi-trait genotype
#'   effects. If \code{NULL}, the trait covariance matrix is estimated from
#'   \code{genotype.effects}.
#' @param desired.gains Optional numeric vector specifying the desired relative
#'   response for each trait. The relative values define the desired direction
#'   of multi-trait genetic improvement.
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{index}}{A data frame containing the selection-index value for
#'     each genotype.}
#'     \item{\code{weights}}{The index coefficients used to construct the
#'     selection index. For a Desired Gains Index, these are derived from the
#'     desired gains and covariance structure; otherwise they are the supplied
#'     coefficients.}
#'   }
#'
#' @references
#' Smith, H. F. (1936).
#' A discriminant function for plant selection.
#' \emph{Annals of Eugenics}, 7(3), 240--250.
#'
#' Hazel, L. N. (1943).
#' The genetic basis for constructing selection indexes.
#' \emph{Genetics}, 28(6), 476--490.
#'
#' Hazel, L. N. and Lush, J. L. (1942).
#' The efficiency of three methods of selection.
#' \emph{Journal of Heredity}, 33(11), 393--399.
#' \doi{10.1093/oxfordjournals.jhered.a105102}
#'
#' Pesek, J. and Baker, R. J. (1969).
#' Comparison of tandem and index selection in the modified pedigree method
#' of breeding self-pollinated species.
#' \emph{Canadian Journal of Plant Science}, 49(6), 773--781.
#' \doi{10.4141/cjps69-132}
#'
#' Werner, C. R., Gardner, K. A. and Tolhurst, D. J. (2024).
#' Reviving the Desired Gains Index: An optimal solution for parent selection
#' in public plant breeding programs.
#' \emph{bioRxiv}, 2024.07.21.603926.
#' \doi{10.1101/2024.07.21.603926}
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
    V.approx <- stats::cov(effects)
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
