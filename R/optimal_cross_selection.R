#' Optimal cross selection via a genetic algorithm (with optional fixed crosses)
#'
#' Uses a genetic algorithm (GA) to choose \code{ncrosses} parent pairs that
#' balance an average criterion (utility, \code{criterion}) and genomic similarity
#' (computed from \code{G.mat}) according to a target trade-off angle.
#' If \code{fixed.crosses} is provided, those crosses are always included and the
#' GA selects the remaining crosses from \code{crosses}.
#'
#' Both \code{crosses} and \code{fixed.crosses} may be provided either as
#' integer indices (1-based, referring to rows/columns of \code{G.mat}) or as
#' character identifiers matching \code{rownames(G.mat)}. I
#'
#' @param crosses Integer or character data.frame or matrix (\eqn{K \times 2}). Candidate
#'   \emph{variable} crosses. Must have exactly two columns; selfings
#'   are not allowed.
#' @param fixed.crosses Integer or character data.frame or matrix (\eqn{F \times 2}) or \code{NULL}.
#'   Crosses that must be included in the final plan. Can have 0 rows. Must use
#'   the same indexing/identifiers as \code{crosses}. Selfings are not allowed.
#' @param ncrosses Integer. Total number of crosses in the final plan
#'   (variable + fixed). Must be \eqn{\ge} \code{nrow(fixed.crosses)}.
#' @param target.angle Numeric, target trade-off angle between utility
#'   and similarity. Smaller angles emphasize \code{criterion}; larger angles
#'   emphasize similarity control, allowed values are in the range of 0 to 90.
#'   Use 0 to prioritize criterion only, and 90 to prioritize diversity (G.mat).
#'   Everything in between is a trade-off.
#' @param criterion Numeric vector (length = \code{nrow(crosses)}). Utility for
#'   each candidate cross in \code{crosses}.
#' @param criterion.fixed Numeric vector (length = \code{nrow(fixed.crosses)}) or
#'   \code{NULL}. Utility for each fixed cross; if \code{fixed.crosses} has 0 rows,
#'   this can be length-0.
#' @param G.mat Numeric square matrix (\eqn{n \times n}). Genomic relationship matrix
#'   used to evaluate similarity/diversity among parents. If \code{crosses} or
#'   \code{fixed.crosses} are character, \code{rownames(G.mat)} must be present.
#' @param params A list of GA parameters:
#'   \describe{
#'     \item{\code{propability.mutate}}{Numeric in \eqn{[0,1]}. Mutation probability.}
#'     \item{\code{n.mutate}}{Integer \eqn{\ge 0}. Mutation magnitude parameter.}
#'     \item{\code{n.select}}{Integer \eqn{\ge 2}. Number selected per generation.}
#'     \item{\code{n.pop}}{Integer \eqn{\ge 2}. Population size.}
#'     \item{\code{max.generation}}{Integer \eqn{\ge 1}. Maximum generations.}
#'     \item{\code{max.iteration}}{Integer \eqn{\ge 1}. Max iterations without improvement.}
#'     \item{\code{angle.penalty}}{Numeric \eqn{\ge 0}. Penalty for deviation from \code{target.angle}.}
#'   }
#' @param return.params Logical. If \code{TRUE}, return additional GA summary metrics.
#' @param nthreads Integer \eqn{\ge 1}. Number of threads to use.
#'
#' @return
#' \describe{
#'   \item{If \code{return.params = FALSE}:}{A data.frame with columns
#'     \code{parent1} and \code{parent2} describing the selected cross plan,
#'     returned in the same representation (integer indices or character IDs)
#'     as supplied in \code{crosses}/\code{fixed.crosses}.}
#'   \item{If \code{return.params = TRUE}:}{A list with elements:
#'     \describe{
#'       \item{\code{crossPlan}}{A data.frame with columns \code{parent1} and
#'         \code{parent2} describing the selected cross plan (mapped back to the
#'         original input representation).}
#'       \item{\code{uMax}}{Maximum attainable average utility when optimizing only utility.}
#'       \item{\code{uMin}}{Minimum attainable average utility under similarity optimization.}
#'       \item{\code{simMax}}{Similarity corresponding to \code{uMax}.}
#'       \item{\code{simMin}}{Similarity corresponding to \code{uMin}.}
#'       \item{\code{uBest}}{Average utility achieved by the optimized cross plan.}
#'       \item{\code{simBest}}{Similarity of the optimized cross plan.}
#'       \item{\code{angleBest}}{Angle (radians) of the optimized solution.}
#'       \item{\code{lenBest}}{Objective-space length of the optimized solution vector.}
#'     }}
#' }
#'
#' @details
#' The returned cross plan is assembled from the selected rows of \code{crosses}
#' and all rows of \code{fixed.crosses}. If fixed crosses also appear among the
#' candidate crosses, duplicates may occur in the final plan (a warning is issued).
#'
#' @export


optimal_cross_selection <- function(crosses,
                                    fixed.crosses = NULL,
                                    ncrosses,
                                    target.angle,
                                    criterion,
                                    criterion.fixed = NULL,
                                    G.mat,
                                    return.params=FALSE,
                                    params = list(),
                                    nthreads = 4L) {

  defaults <-
    list(
      propability.mutate = 0.01,
      n.mutate = 2,
      n.select = 500,
      n.pop = 10000,
      max.generation = 100,
      max.iteration = 100,
      angle.penalty = 0.5
    )
  G <- G.mat
  n.Threads <- nthreads
  params <- modifyList(defaults, params)
  u <- criterion
  u.fixed = criterion.fixed
  ##  Coerce types

  # target.angle is supplied in DEGREES (0..90)
  if (length(target.angle) != 1L || !is.finite(target.angle)) {
    stop("`target.angle` must be a single finite number (in degrees).")
  }
  if (target.angle < 0 || target.angle > 90) {
    stop("`target.angle` must be between 0 and 90 degrees. ",
         "Use 0 to prioritize criterion only, and 90 to prioritize diversity (G.mat).")
  }

  target.angle <- target.angle * pi / 180


  if (!is.matrix(G)) G <- as.matrix(G)
  crosses <- as.matrix(crosses)



  if(!is.null(fixed.crosses)){
    fixed.crosses <- as.matrix(fixed.crosses)

    fixed.crosses <- as.matrix(fixed.crosses)

    if (is.numeric(fixed.crosses) || is.integer(fixed.crosses)) {
      if (any(!is.finite(fixed.crosses))) stop("`fixed.crosses` contains non-finite entries.")
      if (any(fixed.crosses < 1 | fixed.crosses > nrow(G.mat))) {
        stop("Some genotype indices in `fixed.crosses` are outside 1..nrow(G.mat).")
      }
      fixed.crosses2 <- as.data.frame(fixed.crosses)
    } else {
      if (is.null(rownames(G.mat))) stop("Character `fixed.crosses` requires rownames(G.mat).")

      idx <- match(as.vector(fixed.crosses), rownames(G.mat))
      if (anyNA(idx)) stop("Some entries in `fixed.crosses` are not in rownames(G.mat).")

      idx <- matrix(idx, nrow = nrow(fixed.crosses), ncol = ncol(fixed.crosses), byrow = FALSE)
      fixed.crosses2 <- as.data.frame(idx)
    }

    fixed.crosses2 <- as.matrix(fixed.crosses2)

  }


  crosses <- as.matrix(crosses)

  if (is.numeric(crosses) || is.integer(crosses)) {
    if (any(!is.finite(crosses))) stop("`crosses` contains non-finite entries.")
    if (any(crosses < 1 | crosses > nrow(G.mat))) {
      stop("Some genotype indices in `crosses` are outside 1..nrow(G.mat).")
    }
    crosses2 <- as.data.frame(crosses)
  } else {
    if (is.null(rownames(G.mat))) stop("Character `crosses` requires rownames(G.mat).")

    idx <- match(as.vector(crosses), rownames(G.mat))
    if (anyNA(idx)) stop("Some entries in `crosses` are not in rownames(G.mat).")

    idx <- matrix(idx, nrow = nrow(crosses), ncol = ncol(crosses), byrow = FALSE)
    crosses2 <- as.data.frame(idx)
  }

  crosses2 <- as.matrix(crosses2)


  u       <- as.numeric(u)

  if (is.null(fixed.crosses) || nrow(fixed.crosses) == 0) {
    fixed.crosses2 <- crosses2[0, , drop = FALSE]  # 0-row, 2-col
    u.fixed <- u[0]                              # numeric(0)
  }
  fixed.crosses2 <- as.matrix(fixed.crosses2)
  u.fixed <- as.numeric(u.fixed)

  ##  Basic matrix checks
  if (!is.numeric(G)) stop("`G` must be a numeric matrix.")
  if (nrow(G) != ncol(G)) stop("`G` must be square.")
  nInd <- ncol(G)
  if (nInd < 2) stop("`G` must have at least 2 individuals.")

  # Cross matrices must be 2 columns (allow 0 rows for fixed)
  if (ncol(crosses) != 2) stop("`crosses` must have exactly 2 columns.")
  if (!is.null(fixed.crosses) && nrow(fixed.crosses) > 0L && ncol(fixed.crosses) != 2) {
    stop("`fixed.crosses` must have exactly 2 columns (or 0 rows).")
  }


  ##  Index validity & structure
  if (any(!is.finite(crosses2))) stop("`crosses` contains non-finite entries.")
  if (nrow(fixed.crosses2) && any(!is.finite(fixed.crosses2))) stop("`fixed.crosses` contains non-finite entries.")

  if (any(crosses2 < 1 | crosses2 > nInd)) stop("Some indices in `crosses` are outside 1..nrow(G).")
  if (nrow(fixed.crosses2) && any(fixed.crosses2 < 1 | fixed.crosses2 > nInd)) {
    stop("Some indices in `fixed.crosses` are outside 1..nrow(G).")
  }

  if (any(crosses2[,1] == crosses2[,2])) stop("`crosses` contains self-crosses.")
  if (nrow(fixed.crosses2) && any(fixed.crosses2[,1] == fixed.crosses2[,2])) stop("`fixed.crosses` contains self-crosses.")


  # Overlap warning
  if (!is.null(fixed.crosses) && nrow(fixed.crosses) > 0) {

    ov <- paste0(fixed.crosses[,1], "_", fixed.crosses[,2]) %in%
      paste0(crosses[,1], "_", crosses[,2])
    if (any(ov)) {
      warning(sum(ov), " fixed crosses also appear in `crosses`. ",
              "They will still be appended; resulting plan may contain duplicates of those pairs.")
    }
  }

  ##  Vector checks
  if (!all(is.finite(u))) stop("`u` contains non-finite values.")
  if (nrow(crosses) != length(u)) {
    stop("`u` length (", length(u), ") must equal nrow(crosses) (", nrow(crosses), ").")
  }
  if (ncol(crosses) != 2) stop("`crosses` must have exactly 2 columns.")
  if (nrow(fixed.crosses2) > 0L && ncol(fixed.crosses2) != 2) stop("`fixed.crosses` must have exactly 2 columns (or 0 rows).")


  ##  Scalar checks (use n.Threads)
  req_scalar_num <- list(
    ncrosses = ncrosses,
    target.angle = target.angle,
    propability.mutate = params$propability.mutate,
    nmutate = params$n.mutate,
    nselect = params$n.select,
    npop = params$n.pop,
    max.generation = params$max.generation,
    max.iteration = params$max.iteration,
    angle.penalty = params$angle.penalty,
    `n.Threads` = n.Threads
  )
  for (nm in names(req_scalar_num)) {
    val <- req_scalar_num[[nm]]
    if (length(val) != 1L || !is.finite(val)) stop("`", nm, "` must be a single finite value.")
  }
  # integer-ish
  for (nm in c("ncrosses","nmutate","nselect","npop","max.generation","max.iteration","n.Threads")) {
    v <- req_scalar_num[[nm]]
    if (abs(v - round(v)) > .Machine$double.eps^0.5) stop("`", nm, "` must be an integer.")
  }

  ncrosses       <- as.integer(ncrosses)
  nmutate        <- as.integer(params$n.mutate)
  nselect        <- as.integer(params$n.select)
  npop           <- as.integer(params$n.pop)
  max.generation <- as.integer(params$max.generation)
  max.iteration  <- as.integer(params$max.iteration)
  nThreads       <- as.integer(n.Threads)

  if (ncrosses < 1L) stop("`ncrosses` must be >= 1.")
  if (nmutate < 0L) stop("`nmutate` must be >= 0.")
  if (nselect < 2L) stop("`nselect` must be >= 2.")
  if (npop < 2L) stop("`npop` must be >= 2.")
  if (max.generation < 1L) stop("`max.generation` must be >= 1.")
  if (max.iteration  < 1L) stop("`max.iteration` must be >= 1.")
  if (nThreads < 1L) stop("`n.Threads` must be >= 1.")
  if (!is.numeric(target.angle)) stop("`target.angle` must be numeric.")
  if (!is.numeric(params$propability.mutate) || params$propability.mutate < 0 || params$propability.mutate > 1){
    stop("`propability.mutate` must be in [0,1].")}
  if (!is.numeric(params$angle.penalty) || params$angle.penalty < 0){
    stop("`angle.penalty` must be >= 0.")}

  ##  Cross count logic
  fixedCount <- nrow(fixed.crosses2)
  nVar <- ncrosses - fixedCount
  if (nVar < 0L) {
    stop("`ncrosses` (", ncrosses, ") must be >= number of fixed crosses (", fixedCount, ").")
  }
  if (nVar > nrow(crosses)) {
    stop("Not enough variable candidates: need ", nVar, " but only ", nrow(crosses), " in `crosses`.")
  }
  if (nselect > npop) stop("`nselect` must be <= `npop`.")

  ##  Map names to C++
  res <- cpp_optimal_cross_selection(
    Crosses       = crosses2,
    fixedCrosses  = fixed.crosses2,
    nCross        = ncrosses,
    targetAngle   = as.numeric(target.angle),
    u             = u,
    ufixed        = u.fixed,
    G             = G,
    probMut       = as.numeric(params$propability.mutate),
    nMutate       = nmutate,
    nSel          = nselect,
    nPop          = npop,
    maxGen        = max.generation,
    maxRun        = max.iteration,
    anglePenalty  = as.numeric(params$angle.penalty),
    nThreads      = nThreads
  )
  if(return.params){
    res$crossPlan <- as.data.frame(res$crossPlan)[,]
    temp <- rbind(crosses2,fixed.crosses2)
    if (is.null(fixed.crosses) || nrow(fixed.crosses)==0){
      fixed.orig <- crosses[0,,drop=FALSE]
    } else{
        fixed.orig <- fixed.crosses
    }
    temp2 <- rbind(crosses, fixed.orig)


    idx <- match(
      paste(res$crossPlan[,1], res$crossPlan[,2]),
      paste(temp[,1],     temp[,2])
    )
    res$crossPlan <- as.data.frame(temp2[idx,])
    names(res$crossPlan) <-c("parent1","parent2")
    }else{
      temp <- rbind(crosses2,fixed.crosses2)
      if (is.null(fixed.crosses) || nrow(fixed.crosses)==0){
        fixed.orig <- crosses[0,,drop=FALSE]
      } else{
        fixed.orig <- fixed.crosses
      }
      temp2 <- rbind(crosses, fixed.orig)


      idx <- match(
        paste(res$crossPlan[,1], res$crossPlan[,2]),
        paste(temp[,1],     temp[,2])
      )

      res <- as.data.frame(temp2[idx,])
    names(res) <-c("parent1","parent2")
    }

  return(res)
}
