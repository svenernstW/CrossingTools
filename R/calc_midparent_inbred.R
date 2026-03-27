#' Predicted cross means (additive) from specified crosses
#'
#' Compute the expected genomic estimated breeding value (EGEBV) for each proposed
#' cross as the mean of the parental GEBVs (mid-parent value). Optionally, compute
#' a selection index per cross as a weighted sum across traits.
#'
#' This function supports two-way crosses (2 parents) and four-way crosses (4 parents).
#'
#' @param crosses A matrix or data.frame specifying the parents for each proposed cross:
#'   \itemize{
#'     \item two-way crosses: \code{n_crosses x 2}
#'     \item four-way crosses: \code{n_crosses x 4}
#'   }
#'   Entries may be either integer indices referring to rows of \code{marker.mat} (1-based),
#'   or character identifiers matching \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker matrix with genotypes in rows and markers in columns.
#'   The coding must be consistent with the marker effects in \code{effects} (e.g., dosage coding).
#' @param marker.effects Numeric matrix of marker effects with markers in rows and traits in columns.
#'   Must have \code{nrow(effects) == ncol(marker.mat)}.
#' @param weights Optional numeric vector of trait weights of length \code{ncol(effects)}.
#'   If provided, an index value \code{IDX} is computed for each cross as a weighted sum of
#'   predicted EGEBVs across traits.
#' @param nthreads Integer (default 4). Number of threads used by the C++ backend.
#'
#' @return If \code{weights} is \code{NULL}, returns a data.frame with columns:
#'   \itemize{
#'     \item parent identifiers (\code{parent1}, \code{parent2}; optionally \code{parent3}, \code{parent4})
#'     \item predicted cross means per trait: \code{EGEBV1, EGEBV2, ...}
#'   }
#'   If \code{weights} is provided, returns a list with:
#'   \itemize{
#'     \item \code{cross.df}: the data.frame described above
#'     \item \code{index.df}: a data.frame with parent identifiers and a single column \code{IDX}
#'   }
#'
#' @export
#'
calc_midparent_inbred <- function(crosses,  marker.mat, marker.effects,  weights = NULL,
                              nthreads = 4L) {
  n.Threads <- nthreads
  effects <- marker.effects
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


  # ---- Normalize ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)
  if (!is.matrix(effects)) effects <- as.matrix(effects)
  crosses2 <- as.matrix(crosses2)
  calculate.index <- !is.null(weights)

  if (ncol(effects) == 1 | is.null(weights)) { weights <- rep(1,ncol(effects)) }

  # ---- Checks ----
  if (ncol(marker.mat) <= 0L) stop("marker.mat must have markers in columns.")
  if (nrow(effects) != ncol(marker.mat)) stop("effects must have nrow(effects) == ncol(marker.mat).")

  if (calculate.index && is.null(weights)) {
    stop("`weights` is required when calculate.index = TRUE.")
  }
  if (calculate.index && length(weights) != ncol(effects)) {
    stop("`weights` must have length equal to ncol(effects) when calculate.index = TRUE.")
  }
  if (calculate.index) {
    weights <- as.numeric(weights)
    if (any(!is.finite(weights))) stop("`weights` must contain only finite values.")
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


  name_vec <- paste0(rep(c("GEBV"), each = ncol(effects)), seq_len(ncol(effects)))
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
