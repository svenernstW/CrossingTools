#' Calculate expected breeding values for inbred crosses
#'
#' Calculates the expected genomic breeding value of progeny from proposed
#' biparental or four-parent crosses among inbred lines.
#'
#' Under an additive genetic model, the expected breeding value of progeny from
#' a biparental cross equals the average of the parental breeding values, i.e.
#' the mid-parent value (Falconer and Mackay, 2009). For four-parent crosses,
#' the expected breeding value is analogously obtained as the average breeding
#' value of the four contributing parents.
#'
#' Breeding values are calculated from centred marker genotypes,
#' \eqn{M - 2p}, and average allele-substitution effects (\eqn{\alpha}), where
#' \eqn{p} contains the allele frequencies of the reference population.
#' The reference allele frequencies can be supplied by the user. If they are
#' not supplied, they are calculated from \code{marker.mat}. Thus, \code{p}
#' defines the reference population relative to which breeding values are
#' expressed.
#'
#' Multiple traits can be evaluated simultaneously. Optional trait weights can
#' be used to calculate a weighted breeding-value index for each cross.
#'
#' @param crosses Matrix or data frame with two columns for biparental crosses
#'   or four columns for four-parent crosses. Parent identifiers may be row
#'   indices of \code{marker.mat} or character identifiers matching
#'   \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker dosage matrix with individuals in rows and
#'   markers in columns, coded 0, 1, and 2 for the counted allele.
#' @param marker.effects Numeric matrix of average allele-substitution effects
#'   (\eqn{\alpha}) with markers in rows and traits in columns. Effects should
#'   be parameterised relative to the reference population defined by
#'   \code{p}. If \code{p = NULL}, the reference allele frequencies are derived
#'   from \code{marker.mat}. Its number of rows must equal
#'   \code{ncol(marker.mat)}.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, a weighted genomic breeding-value index is calculated for each
#'   cross.
#' @param p Optional numeric vector of reference allele frequencies, with one
#'   value per marker. These frequencies define the reference population used
#'   to centre marker genotypes for breeding-value calculations. If
#'   \code{NULL}, allele frequencies are calculated from \code{marker.mat}.
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
#' @references
#' Falconer, D. S. and Mackay, T. (2009).
#' \emph{Introduction to Quantitative Genetics}. 4th ed.
#' Pearson, Prentice Hall, Harlow.
#'
#' @export

calc_midparent_inbred <- function(
    crosses,
    marker.mat,
    marker.effects,
    weights = NULL,
    p = NULL,
    nthreads = 4L
) {
  n.Threads <- nthreads

  effects <- as.matrix(marker.effects)

  traits <- colnames(effects)
  if (is.null(traits)) {
    traits <- paste0("trait", seq_len(ncol(effects)))
  }
  if(!ncol(crosses) %in% c(2,4)){stop("`crosses` must have 2 columns for biparental crosses or 4 columns for four-parent crosses.")}
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

  ###############################################################################
  # Reference allele frequencies
  ###############################################################################

  if (is.null(p)) {

    # Default: the supplied marker population defines the reference population
    p <- colMeans(marker.mat) / 2

  } else {

    if (!is.numeric(p)) {
      stop("`p` must be a numeric vector of reference allele frequencies.")
    }

    p <- as.numeric(p)

    if (length(p) != ncol(marker.mat)) {
      stop(
        "`p` must contain one allele frequency per marker: ",
        "length(p) = ", length(p),
        ", ncol(marker.mat) = ", ncol(marker.mat), "."
      )
    }
  }

  if (any(!is.finite(p))) {
    stop("`p` must contain only finite values.")
  }

  if (any(p < 0 | p > 1)) {
    stop("All values in `p` must be between 0 and 1.")
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
        Crosses   = crosses2,
        M         = marker.mat,
        U         = effects,
        weights   = weights,
        p         = p,
        calcindex = calculate.index,
        nThreads  = nThreads
      )

        }


  if(cross.type=="4W"){

    temp <- cpp_calculate_expectation_A4W(
      Crosses   = crosses2,
      M         = marker.mat,
      U         = effects,
      weights   = weights,
      p         = p,
      calcindex = calculate.index,
      nThreads  = nThreads
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
