#' Calculate predicted additive means for specified crosses
#'
#' Calculates the expected genomic breeding value of proposed two-way or
#' four-way crosses as the mean of the parental genomic breeding values.
#'
#' Multiple traits can be evaluated simultaneously. Optional trait weights can
#' be used to calculate a weighted genomic breeding value index for each cross.
#'
#' @param crosses Matrix or data frame with two columns for two-way crosses or
#'   four columns for four-way crosses. Parent identifiers may be row indices of
#'   \code{marker.mat} or character identifiers matching
#'   \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker dosage matrix with individuals in rows and
#'   markers in columns, coded 0, 1, and 2 for the counted allele.
#' @param marker.effects Numeric matrix of average allele-substitution effects
#'   (\eqn{\alpha}) with markers in rows and traits in columns. Effects must be
#'   parameterised relative to the reference population represented by
#'   \code{marker.mat}. Its number of rows must equal
#'   \code{ncol(marker.mat)}.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, a weighted genomic breeding value index is calculated for each
#'   cross.
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return If \code{weights = NULL}, a data frame containing the parental
#'   identifiers followed by one \code{GEBV.<trait>} column per trait.
#'
#'   If \code{weights} is supplied, a list containing:
#'   \describe{
#'     \item{\code{cross.df}}{The parental identifiers and trait-specific
#'     genomic breeding values.}
#'     \item{\code{index.df}}{The parental identifiers and the weighted index
#'     \code{GEBV.IDX}.}
#'   }
#'
#' @export

calc_midparent_inbred <- function(crosses,  marker.mat, marker.effects,  weights = NULL,
                              nthreads = 4L) {
  n.Threads <- nthreads

  effects <- as.matrix(marker.effects)

  traits <- colnames(effects)
  if (is.null(traits)) {
    traits <- paste0("trait", seq_len(ncol(effects)))
  }
  if(!ncol(crosses) %in% c(2,4)){stop("ncol(crosses) needs to be 2 for two way crosses or 4 for three or four way crosses")}
  crosses_in <- crosses

  crosses <- as.matrix(crosses)

  if (is.numeric(crosses) || is.integer(crosses)) {
    if (any(!is.finite(crosses))) stop("`crosses` contains non-finite entries.")
    if (any(crosses < 1 | crosses > nrow(marker.mat))) {
      stop("Some genotype indices in `crosses` are outside 1..nrow(marker.mat).")
    }
    crosses2 <- as.data.frame(crosses)
  } else {
    if (is.null(rownames(marker.mat))) stop("Character `crosses` requires rownames(marker.mat).")

    idx <- match(as.vector(crosses), rownames(marker.mat))
    if (anyNA(idx)) stop("Some entries in `crosses` are not in rownames(marker.mat).")

    idx <- matrix(idx, nrow = nrow(crosses), ncol = ncol(crosses), byrow = FALSE)
    crosses2 <- as.data.frame(idx)
  }



  if (ncol(crosses) == 2){
    cross.type <- "2W"
  }  else {
    cross.type <- "4W"}


  if (cross.type == "2W") {
    names(crosses2) <- c("Parent1","Parent2")
  } else {
    names(crosses2) <- c("Parent1","Parent2","Parent3","Parent4")
  }

  crosses2 <- as.matrix(crosses2)
  storage.mode(crosses2) <- "integer"
  if (anyNA(crosses2)) stop("Internal error: `crosses2` contains NA after conversion.")

  if (!is.numeric(marker.mat)) {
    stop("`marker.mat` must be numeric.")
  }

  if (any(!is.finite(marker.mat))) {
    stop("`marker.mat` must contain only finite values.")
  }

  if (any(!marker.mat %in% c(0, 1, 2))) {
    stop("`marker.mat` must contain genotype dosages 0, 1, or 2.")
  }

  if (!is.numeric(effects)) {
    stop("`marker.effects` must be numeric.")
  }

  if (any(!is.finite(effects))) {
    stop("`marker.effects` must contain only finite values.")
  }

  # ---- Normalize ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  if (!is.matrix(effects)) effects <- as.matrix(effects)
  crosses2 <- as.matrix(crosses2)
  calculate.index <- !is.null(weights)

  # ---- Checks ----
  if (ncol(marker.mat) <= 0L) {
    stop("`marker.mat` must have markers in columns.")
  }

  if (nrow(effects) != ncol(marker.mat)) {
    stop(
      "`marker.effects` must have one row per marker: ",
      "nrow(marker.effects) = ", nrow(effects),
      ", ncol(marker.mat) = ", ncol(marker.mat), "."
    )
  }

  if (is.null(weights)) {
    # Dummy value required by C++; not used when calcindex = FALSE
    weights <- rep(1, ncol(effects))
  } else {
    weights <- as.numeric(weights)

    if (length(weights) != ncol(effects)) {
      stop("`weights` must have length equal to ncol(marker.effects).")
    }

    if (any(!is.finite(weights))) {
      stop("`weights` must contain only finite values.")
    }
  }






  if (length(n.Threads) != 1L || !is.finite(n.Threads) || n.Threads < 1 || n.Threads != as.integer(n.Threads)) {
    stop("`n.Threads` must be a positive integer.")
  }
  nThreads <- as.integer(n.Threads)

    if(cross.type=="2W"){

        temp <- cpp_calculate_expectation_A(
          Crosses    = crosses2,
           M = marker.mat,
          U    = effects,
           weights    = weights,
          calcindex  = calculate.index,
          nThreads   = nThreads
        )
        }


  if(cross.type=="4W"){

    temp <- cpp_calculate_expectation_A4W(
      Crosses    = crosses2,
      M = marker.mat,
      U    = effects,
      weights    = weights,
      calcindex  = calculate.index,
      nThreads   = nThreads
    )
  }


  name_vec <- paste0("GEBV.",traits)
  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- if (ncol(crosses) == 2) c("parent1","parent2") else c("parent1","parent2","parent3","parent4")


  if (calculate.index) {
    temp1 <- as.data.frame(temp)[, 1:ncol(effects), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    temp2 <- as.data.frame(temp)[, ((ncol(effects)) + 1):ncol(as.data.frame(temp)), drop = FALSE]
    names(temp2) <- c("GEBV.IDX")
    temp2 <- cbind(crosses_df, temp2)
    return(list(cross.df = temp1, index.df = temp2))
  } else {
    temp1 <- as.data.frame(temp)[, 1:( ncol(effects)), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)
    return(temp1)
  }
}
