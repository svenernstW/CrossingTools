#' Summarize a multi-trait selection index
#'
#' Computes summary statistics for a linear selection index given a trait covariance
#' matrix and either (i) desired gains or (ii) fixed weights (Smith-Hazel).
#'
#' For a trait covariance matrix \eqn{G} and index weights \eqn{b}, the index is
#' \eqn{I = b^{\top} T}. The function returns the index variance
#' \eqn{\sigma_I^2 = b^{\top} G b}, the expected standardized response in the index
#' \eqn{i \sigma_I} (where \code{i} is the standardized selection differential),
#' the expected response in each trait \eqn{i \, G b / \sigma_I}, and the correlation
#' between each trait and the index.
#'
#' If \code{gains} are provided, weights are computed as \eqn{b = G^{-1} d}, where
#' \eqn{d} is the desired gains vector.
#'
#' @param var.mat Numeric trait covariance matrix (\eqn{n_{trait} \times n_{trait}}).
#' @param gains Optional numeric vector of desired gains (length \eqn{n_{trait}}).
#'   If supplied, index weights are computed as \code{solve(var.mat) \%*\% gains}.
#' @param weights Optional numeric vector of index weights (length \eqn{n_{trait}}).
#'   Supply this to summarize an index with fixed weights.Using cov(EBVs) (i.e. cov(BLUPs))
#'   as var.mat provides an empirical implementation of the Henderson form of the Smith–Hazel index.
#' @param intensity Numeric scalar. Standardized selection differential \eqn{i} used to scale expected responses
#'   (default 1).
#' @param plot Logical. If \code{TRUE}, produce a bar plot of expected trait gains
#'   (default \code{TRUE}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{overall.df}: data.frame with \code{index.var} (\eqn{\sigma_I^2}) and \code{gain.index} (\eqn{i\sigma_I})
#'   \item \code{traits.df}: data.frame with per-trait variance, correlation with the index, expected gain, and weight
#' }
#'
#' @export



measure_traits <- function(var.mat=NULL, gains = NULL, weights = NULL, intensity=1,plot=TRUE) {
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
    if(ncol(var.mat)!= length(gains)){
      stop("gains needs to have the same length as ncol(var.mat)")
    }

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
    if(ncol(var.mat)!= length(weights)){
      stop("weights needs to have the same length as ncol(var.mat)")
    }
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
