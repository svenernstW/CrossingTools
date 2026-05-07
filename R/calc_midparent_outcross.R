#' Predicted outcross means from specified two-way crosses (EGEBV and ETGV)
#'
#' Compute predicted cross means for two-way (bi-parental) crosses using marker effects.
#' For each cross and each trait, the function returns:
#' \itemize{
#'   \item \code{GEBV}: the mid-parent mean of parental additive genomic values.
#'   \item \code{TGV}: the predicted mean total genotypic value of the F1, combining additive
#'         and (optional) dominance contributions (Werner et al).
#' }
#' If \code{weights} are provided, two index values are also returned: a weighted index based
#' on the \code{EGEBV} trait means and a weighted index based on the \code{ETGV} trait means.
#'
#' @param crosses A matrix or data.frame (\code{n_crosses x 2}) specifying the parents for each
#'   proposed cross. Entries may be either integer indices referring to rows of \code{marker.mat},
#'   or character identifiers matching \code{rownames(marker.mat)}.
#' @param marker.mat Numeric marker dosage matrix with genotypes in rows and markers in columns.
#'   The coding must be consistent with the marker effects provided in \code{effects.A} and
#'   \code{effects.D}.
#' @param marker.effects.A Numeric matrix of additive marker effects (markers in rows, traits in columns).
#'   Must satisfy \code{nrow(marker.effects.A) == ncol(marker.mat)}.
#' @param marker.effects.D Optional numeric matrix of dominance marker effects (markers in rows, traits in columns).
#'   Must have the same dimensions as \code{marker.effects.A}. If \code{NULL}, dominance is ignored (treated as zero).
#' @param weights Optional numeric vector of trait weights of length \code{ncol(marker.effects.A)}. If provided,
#'   index values \code{GEBV.IDX} and \code{TGV.IDX} are computed per cross as weighted sums across traits.
#' @param nthreads Integer (default 4). Number of threads used by the C++ backend.
#'
#' @return If \code{weights} is \code{NULL}, a data.frame with columns:
#'   \itemize{
#'     \item \code{parent1}, \code{parent2}
#'     \item predicted trait means: \code{EGEBV1..EGEBVt}, \code{ETGV1..ETGVt}
#'   }
#'   If \code{weights} is provided, a list with:
#'   \itemize{
#'     \item \code{cross.df}: the data.frame described above
#'     \item \code{index.df}: a data.frame with \code{parent1}, \code{parent2}, \code{GEBV.IDX}, \code{TGV.IDX}
#'   }
#'
#' @export

#'
calc_midparent_outcross <- function(crosses,  marker.mat, marker.effects.A, marker.effects.D=NULL,  weights = NULL,
                                   nthreads = 4L) {

  nThreads = as.integer(nthreads)
  hap.mat1 <- NULL
  hap.mat2 <- NULL
  effects.A <- as.matrix(marker.effects.A)
  effects.D <- as.matrix(marker.effects.D)

  if (is.null(marker.mat)) {
    stop("marker.mat must be provided.")
  }


  if(is.null(marker.mat) & !is.null(hap.mat1) & !is.null(hap.mat2)){
    print("marker.mat not provided, calculating it from hap.mat1 and hap.mat2")
    marker.mat <- hap.mat1+hap.mat2
  }

  # ---- Normalize ----
  if (!is.matrix(marker.mat)) marker.mat <- as.matrix(marker.mat)

  calculate.index <- !is.null(weights)

  if (is.null(effects.D) ) {
    warning("No effects.D given, ignoring dominane effects")
    effects.D <- matrix(0,ncol=ncol(effects.A), nrow=nrow(effects.A))
  }
  if (!is.matrix(effects.A)) effects.A <- as.matrix(effects.A)
  if (!is.matrix(effects.D)) effects.D <- as.matrix(effects.D)

  if (ncol(effects.A) == 1 | is.null(weights)) { weights <- rep(1,ncol(effects.A)) }
  if (!all(dim(effects.D) == dim(effects.A))) {
    stop("effects.A and effects.D need to have the same dimensions.")
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


  name_vec <- paste0(rep(c("GEBV","TGV"), each = ncol(effects.A)), seq_len(ncol(effects.A)))
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
