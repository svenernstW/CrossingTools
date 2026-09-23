#' Predict response to multi-trait index selection
#'
#' Calculates the expected response to multi-trait selection for either a
#' desired-gains index or an index with user-defined trait weights.
#'
#' Linear selection indices combine information from multiple traits into a
#' single selection criterion and provide a framework for predicting correlated
#' responses to selection (Smith, 1936; Hazel, 1943). Index selection can
#' improve multiple traits simultaneously and has long been applied in plant
#' and animal breeding (Hazel and Lush, 1942; Pesek and Baker, 1969).
#'
#' For a desired-gains index, index weights are derived from a vector of desired
#' responses \eqn{g} as
#' \deqn{b = \Sigma^{-1} g,}
#' where \eqn{\Sigma} is the trait covariance matrix. The desired-gains
#' formulation specifies the relative direction of genetic improvement rather
#' than fixing the absolute magnitude of response. Consequently, the predicted
#' responses are proportional to the supplied desired gains, while their
#' overall magnitude depends on the covariance structure and the selection
#' intensity. This formulation allows breeding objectives to be expressed
#' directly as desired relative changes among traits without requiring explicit
#' economic weights (Werner et al., 2024).
#'
#' Alternatively, fixed index weights can be supplied directly. This corresponds
#' to the classical linear selection-index framework of Smith (1936) and
#' Hazel (1943), in which multiple traits are combined into a single selection
#' criterion. Here, the supplied weights are treated directly as the index
#' coefficients \eqn{b}. For a vector of index coefficients \eqn{b}, the index
#' variance is
#' \deqn{\sigma_I^2 = b^\top \Sigma b,}
#' and the expected response in trait \eqn{k} is
#' \deqn{\Delta G_k =
#' i \frac{(\Sigma b)_k}{\sqrt{b^\top \Sigma b}},}
#' where \eqn{i} is the selection intensity. The function additionally returns
#' the correlation between each trait and the selection index.
#'
#' Exactly one of \code{desired.gains} or \code{weights} must be supplied.
#'
#' If \code{var.mat} is a block-structured covariance matrix containing one
#' trait covariance block per genotype, the diagonal trait covariance blocks
#' are averaged before constructing the index. Covariances between different
#' genotypes are therefore not used in the response calculation.
#'
#' @param var.mat Numeric square covariance matrix. It may either be an
#'   \eqn{n_{\mathrm{trait}} \times n_{\mathrm{trait}}} trait covariance
#'   matrix or a block-structured covariance matrix containing one
#'   \eqn{n_{\mathrm{trait}} \times n_{\mathrm{trait}}} diagonal block per
#'   genotype. In the latter case, the diagonal blocks are averaged to obtain
#'   the trait covariance matrix used for index calculations.
#' @param desired.gains Optional numeric vector specifying the desired relative
#'   response for each trait. Index weights are calculated from the desired
#'   response vector and the trait covariance matrix as
#'   \eqn{b = \Sigma^{-1}g}.
#' @param weights Optional numeric vector specifying the fixed Smith--Hazel
#'   index coefficients for each trait.
#' @param intensity Numeric scalar giving the standardised selection intensity
#'   used to scale the expected responses. The default is 1.
#' @param plot Logical. If \code{TRUE}, plot the expected response for each
#'   trait. The default is \code{TRUE}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{overall.df}}{A data frame containing the index variance
#'     \code{index.var} and expected response in the index
#'     \code{gain.index}.}
#'     \item{\code{trait.df}}{A data frame containing, for each trait, its
#'     variance, correlation with the index, expected response, and index
#'     weight.}
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


predict_response <- function(var.mat=NULL, desired.gains = NULL, weights = NULL, intensity=1,plot=TRUE) {
  gains <- desired.gains
  if(nrow(var.mat) != ncol(var.mat)){
    stop("`var.mat` must be square.")
  }



  if(!is.null(gains) & !is.null(weights) ){
    stop("Provide either weights or gains, function can only handle one at a time!")
  }

  if(is.null(gains) & is.null(weights) ){
    stop("Provide either weights or gains!")
  }




  if(!is.null(gains)){
    ntraits <- length(gains)
    check <- ncol(var.mat)/ntraits

    if(check!= round(check)){
      stop("for desired.gains var.mat needs to be ntrait x ntrait or ntrait*nindividual x ntrait*nindividual")
    }

    temp.var <- matrix(0,nrow = ntraits,ncol = ntraits)

    for(i in 1:check){
      idx <- ((i-1)*ntraits + 1):(i*ntraits)
      temp.var <- temp.var + var.mat[idx, idx]
    }

    temp.var <- temp.var / check

    var.mat <- temp.var

    if(is.null(row.names(var.mat))) row.names(var.mat) <- 1: nrow(var.mat)
    dg_weights <- as.vector(solve(var.mat) %*% gains)

    var_index <- as.vector(t(dg_weights) %*% var.mat %*% dg_weights)

    gain <- intensity * sqrt(var_index)

    overall.df <- data.frame(index.var=var_index,gain.index=gain)

    var_traits <- diag(var.mat)

    gain_traits <- intensity * (var.mat %*% dg_weights)/sqrt(var_index)

    cor_traits <- as.vector(var.mat %*% dg_weights)/(sqrt(var_index) * sqrt(var_traits))

    trait.df <- data.frame(trait=factor(row.names(var.mat)),var = var_traits,cor.index=cor_traits,gain=gain_traits,weight=dg_weights)

  }

  if(!is.null(weights)){
    ntraits <- length(weights)
    check <- ncol(var.mat)/ntraits

    if(check!= round(check)){
      stop("for weights var.mat needs to be ntrait x ntrait or ntrait*nindividual x ntrait*nindividual")
    }

    temp.var <- matrix(0,nrow = ntraits,ncol = ntraits)

    for(i in 1:check){
      idx <- ((i-1)*ntraits + 1):(i*ntraits)
      temp.var <- temp.var + var.mat[idx, idx]
    }

    temp.var <- temp.var / check

    var.mat <- temp.var
    if(is.null(row.names(var.mat))) row.names(var.mat) <- 1: nrow(var.mat)
    var_index <- as.vector(t(weights) %*% var.mat %*% weights)

    gain <- intensity * sqrt(var_index)

    overall.df <- data.frame(index.var=var_index,gain.index=gain)

    var_traits <- diag(var.mat)

    gain_traits <- intensity * (var.mat %*% weights)/sqrt(var_index)

    cor_traits <- as.vector(var.mat %*% weights)/(sqrt(var_index) * sqrt(var_traits))

    trait.df <- data.frame(trait=factor(row.names(var.mat)),var = var_traits,cor.index=cor_traits,gain=gain_traits,weight=weights)

  }
  if(plot){

    p <- ggplot2::ggplot(trait.df, ggplot2::aes(x = gain, y = trait)) +
      ggplot2::geom_col(fill = "white", color = "black", linewidth = 0.6) +
      ggplot2::labs(x = "Gain", y = "Trait") +
      ggplot2::theme_grey(base_size = 10) +
      ggplot2::theme(
        legend.position   = "bottom",
        legend.key        = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        panel.background  = ggplot2::element_rect(
          colour = "black", fill = "grey93", linewidth = 1.1
        ),
        axis.title        = ggplot2::element_text(size = 11),
        axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
        axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
        axis.ticks        = ggplot2::element_line(),
        axis.text         = ggplot2::element_text(size = 10)
      )

    print(p)

  }

  return(list(overall.df = overall.df,trait.df=trait.df))



}
