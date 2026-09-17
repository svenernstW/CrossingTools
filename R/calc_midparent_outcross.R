#' Calculate predicted means for two-way outcrosses
#'
#' Calculates the expected genomic breeding value and expected total genetic
#' value of proposed two-way crosses from additive and dominance marker effects.
#'
#' The expected genomic breeding value is the mid-parent additive value. The
#' expected total genetic value additionally includes the expected dominance
#' contribution of the F1. If dominance effects are not supplied, they are
#' assumed to be zero.
#'
#' Multiple traits can be evaluated simultaneously. Optional trait weights can
#' be used to calculate additive and total-genetic-value indices.
#'
#' @param crosses Matrix or data frame with two columns specifying the parents
#'   of each proposed cross. Parent identifiers may be row indices of
#'   \code{marker.mat} or character identifiers matching
#'   \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker dosage matrix with individuals in rows and
#'   markers in columns.
#' @param marker.effects.A Numeric matrix of average allele-substitution effects
#'   (\eqn{\alpha}) with markers in rows and traits in columns. Effects must be
#'   parameterised relative to the reference population represented by
#'   \code{hap.mat1} and \code{hap.mat2}. Its number of rows must equal
#'   \code{ncol(hap.mat1)}.
#' @param marker.effects.D Numeric matrix of dominance effects with markers
#'   in rows and traits in columns, used with the statistical dominance-deviation
#'   parameterisation. It must have the same dimensions as
#'   \code{marker.effects.A}.
#' @param weights Optional numeric vector with one weight per trait. When
#'   supplied, weighted indices are calculated from the trait-specific GEBVs
#'   and TGVs.
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return If \code{weights = NULL}, a data frame containing the parental
#'   identifiers followed by \code{GEBV.<trait>} and \code{TGV.<trait>} for
#'   each trait.
#'
#'   If \code{weights} is supplied, a list containing:
#'   \describe{
#'     \item{\code{cross.df}}{The parental identifiers and trait-specific GEBVs
#'     and TGVs.}
#'     \item{\code{index.df}}{The parental identifiers and the weighted indices
#'     \code{GEBV.IDX} and \code{TGV.IDX}.}
#'   }
#'
#' @export


calc_midparent_outcross <- function(crosses,  marker.mat, marker.effects.A, marker.effects.D=NULL,  weights = NULL,
                                   nthreads = 4L) {

  effects.A <- as.matrix(marker.effects.A)

  traits <- colnames(effects.A)

  if (is.null(traits)) {
    traits <- paste0("trait", seq_len(ncol(effects.A)))
  }

  if (length(nthreads) != 1L ||
      !is.numeric(nthreads) ||
      !is.finite(nthreads) ||
      nthreads < 1 ||
      nthreads != as.integer(nthreads)) {
    stop("`nthreads` must be a positive integer.")
  }

  nThreads <- as.integer(nthreads)
  nThreads = as.integer(nthreads)
  hap.mat1 <- NULL
  hap.mat2 <- NULL
  effects.A <- as.matrix(marker.effects.A)

  if (is.null(marker.effects.D)) {
    warning("No marker.effects.D supplied; dominance effects are assumed to be zero.",
            call. = FALSE)

    effects.D <- matrix(
      0,
      nrow = nrow(effects.A),
      ncol = ncol(effects.A)
    )

    colnames(effects.D) <- colnames(effects.A)

  } else {
    effects.D <- as.matrix(marker.effects.D)
  }
  if (is.null(marker.mat)) {
    stop("marker.mat must be provided.")
  }



  # ---- Normalize ----
  marker.mat <- as.matrix(marker.mat)

  if (!is.numeric(marker.mat)) {
    stop("`marker.mat` must be numeric.")
  }

  if (any(!is.finite(marker.mat))) {
    stop("`marker.mat` must contain only finite values.")
  }

  if (!all(dim(effects.D) == dim(effects.A))) {
    stop("`marker.effects.A` and `marker.effects.D` must have the same dimensions.")
  }

  # ---- Checks ----
  if (ncol(marker.mat) <= 0L) {
    stop("`marker.mat` must have markers in columns.")
  }

  if (nrow(effects.A) != ncol(marker.mat)) {
    stop(
      "`marker.effects.A` must have one row per marker: ",
      "nrow(marker.effects.A) = ", nrow(effects.A),
      ", ncol(marker.mat) = ", ncol(marker.mat), "."
    )
  }

  calculate.index <- !is.null(weights)

  if (is.null(weights)) {

    # Dummy weights required by the C++ interface.
    # They are not used because calcindex = FALSE.
    weights <- rep(1, ncol(effects.A))

  } else {

    weights <- as.numeric(weights)

    if (length(weights) != ncol(effects.A)) {
      stop(
        "`weights` must have length equal to the number of traits in ",
        "`marker.effects.A`."
      )
    }

    if (any(!is.finite(weights))) {
      stop("`weights` must contain only finite values.")
    }
  }

  # ---- Checks ----
  if (ncol(marker.mat) <= 0L) stop("marker.mat must have markers in columns.")
  if (nrow(effects.A) != ncol(marker.mat)) stop("effects must have nrow(effects) == ncol(marker.mat).")

  if (calculate.index && is.null(weights)) {
    stop("`weights` is required when calculate.index = TRUE.")
  }
  if (calculate.index && length(weights) != ncol(effects.A)) {
    stop("`weights` must have length equal to ncol(effects) when calculate.index = TRUE.")
  }
  if (calculate.index) {
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) stop("`weights` must contain only finite values.")
  }

  if(!ncol(crosses) %in% c(2)){stop("ncol(crosses) needs to be 2")}
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




  crosses2 <- as.matrix(crosses2)
  storage.mode(crosses2) <- "integer"
  if (anyNA(crosses2)) stop("Internal error: `crosses2` contains NA after conversion.")




  Hap1 <- marker.mat/2
  Hap2 <- marker.mat/2



    temp <- cpp_calculate_expectation_AD(
      Crosses    = crosses2,
      Hap1 = Hap1,
      Hap2 = Hap2,
      U    = effects.A,
      D    = effects.D,
      weights    = weights,
      calcindex  = calculate.index,
      nThreads   = nThreads
)


  name_vec <- paste0(rep(c("GEBV.","TGV."), each = ncol(effects.A)),traits)
  crosses_df <- as.data.frame(crosses_in, stringsAsFactors = FALSE)
  names(crosses_df) <- if (ncol(crosses) == 2) c("parent1","parent2") else c("parent1","parent2","parent3","parent4")


  if (calculate.index) {
    temp1 <- as.data.frame(temp)[, 1:(2*ncol(effects.A)), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)

    temp2 <- as.data.frame(temp)[, ((2*ncol(effects.A)) + 1):ncol(as.data.frame(temp)), drop = FALSE]
    names(temp2) <- c("GEBV.IDX","TGV.IDX")
    temp2 <- cbind(crosses_df, temp2)
    return(list(cross.df = temp1, index.df = temp2))
  } else {
    temp1 <- as.data.frame(temp)[, 1:(2* ncol(effects.A)), drop = FALSE]
    names(temp1) <- name_vec
    temp1 <- cbind(crosses_df, temp1)
    return(temp1)
  }
}
