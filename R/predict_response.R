#' Predict the response to multi-trait index selection
#'
#' Summarizes a linear multi-trait selection index from a trait covariance
#' matrix and either desired gains or fixed trait weights.
#'
#' For a desired-gains index, trait weights are derived from the supplied
#' desired gains. For an index with fixed weights, the supplied weights are
#' used directly. Exactly one of \code{desired.gains} or \code{weights} must be
#' provided.
#'
#' The function returns the index variance, the expected response in the index,
#' and the expected response and correlation with the index for each trait.
#'
#' @param var.mat Numeric covariance matrix. It may be either a trait covariance
#'   matrix or a block-structured covariance matrix containing one trait
#'   covariance block per genotype. In the latter case, the trait covariance
#'   matrix is obtained by averaging the diagonal blocks.
#' @param desired.gains Optional numeric vector specifying the desired gain for
#'   each trait.
#' @param weights Optional numeric vector specifying the index weight for each
#'   trait.
#' @param intensity Numeric scalar giving the standardized selection intensity
#'   used to scale the expected responses. The default is 1.
#' @param plot Logical. If \code{TRUE}, plot the expected response for each
#'   trait. The default is \code{TRUE}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{overall.df}}{A data frame containing the index variance
#'     \code{index.var} and expected index response \code{gain.index}.}
#'     \item{\code{trait.df}}{A data frame containing, for each trait, its
#'     variance, correlation with the index, expected response, and index
#'     weight.}
#'   }
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
