#' Optimal cross selection via a genetic algorithm (with optional fixed crosses)
#'
#' Uses a simple genetic algorithm (GA) to choose \code{ncrosses} parent pairs
#' that balance an average criterion (utility \code{u}) and genomic similarity
#' according to a target trade-off angle. Fixed crosses (if any) are forced
#' into the final solution, and the GA fills the remaining slots from
#' \code{crosses}.
#'
#' @param crosses integer matrix (K x 2). Candidate \emph{variable} crosses (1-based indices
#'   referring to rows/columns of \code{G}). Must have exactly two columns
#'   (P1, P2); selfings are not allowed.
#' @param fixed.crosses integer matrix (F x 2) or \code{NULL}. Crosses that must be
#'   included in the final plan. Can have 0 rows. Must use the same indexing as
#'   \code{crosses}. Selfings are not allowed.
#' @param ncrosses integer. Total number of crosses in the final plan
#'   (variable + fixed). Must be \eqn{\ge} number of rows in \code{fixed.crosses}.
#' @param target.angle numeric (radians). Target trade-off angle between the
#'   average criterion and the similarity objective. Smaller angles emphasize
#'   utility \code{u}; larger angles emphasize similarity control.
#' @param u numeric vector (length = nrow(crosses)). Utility (average criterion)
#'   for each candidate cross in \code{crosses}.
#' @param u.fixed numeric vector (length = nrow(fixed.crosses)) or \code{NULL}.
#'   Utility for each fixed cross; if \code{fixed.crosses} has 0 rows, this can
#'   be length-0.
#' @param G numeric square matrix (nInd x nInd). Genomic relationship matrix
#'   used to evaluate similarity/diversity among parents; indices in
#'   \code{crosses}/\code{fixed.crosses} must lie in \code{1..nrow(G)}.
#' @param params A list of parameters for the Genetic Algorithm optimization:
#'   \describe{
#'     \item{\code{propability.mutate}}{numeric in \eqn{[0,1]}. Probability that a candidate solution mutates in the GA.}
#'     \item{\code{n.mutate}}{integer \eqn{\ge 0}. Number of point mutations to apply when a mutation occurs.}
#'     \item{\code{n.select}}{integer \eqn{\ge 2}. Number of candidate solutions selected per generation (selection pool size).}
#'     \item{\code{n.pop}}{integer \eqn{\ge 2}. Population size (number of candidate solutions) per generation.}
#'     \item{\code{max.generation}}{integer \eqn{\ge 1}. Maximum number of GA generations.}
#'     \item{\code{max.iteration}}{integer \eqn{\ge 1}. Maximum number of iterations without improvement.}
#'     \item{\code{angle.penalty}}{numeric \eqn{\ge 0}. Penalty applied to solutions that deviate from \code{target.angle}.}
#'   }
#' @param return.params logical. If \code{TRUE}, return additional information from the GA.
#' @param n.Threads integer \eqn{\ge 1}. Number of OpenMP threads to use.
#'
#' @return
#' \describe{
#'   \item{If \code{return.params = FALSE}:}{A data.frame with columns \code{parent1} and \code{parent2} describing the selected cross plan.}
#'   \item{If \code{return.params = TRUE}:}{A list with elements:
#'     \describe{
#'       \item{\code{crossPlan}}{A data.frame with columns \code{parent1} and \code{parent2} describing the selected cross plan.}
#'       \item{\code{uMax}}{Numeric scalar: maximum attainable average criterion when optimizing only \code{u}.}
#'       \item{\code{uMin}}{Numeric scalar: minimum attainable average criterion when optimizing only similarity.}
#'       \item{\code{simMax}}{Numeric scalar: similarity corresponding to \code{uMax}.}
#'       \item{\code{simMin}}{Numeric scalar: similarity corresponding to \code{uMin}.}
#'       \item{\code{uBest}}{Numeric scalar: average criterion achieved by the optimized cross plan.}
#'       \item{\code{simBest}}{Numeric scalar: similarity of the optimized cross plan.}
#'       \item{\code{angleBest}}{Numeric scalar (radians): angle of the optimized solution relative to the target trade-off.}
#'       \item{\code{lenBest}}{Numeric scalar: length (magnitude) of the optimized solution vector in the objective space.}
#'     }}
#' }
#'
#' @export

optimal_cross_selection <- function(crosses,
                                    fixed.crosses = NULL,
                                    ncrosses,
                                    target.angle,
                                    u,
                                    u.fixed = NULL,
                                    G,
                                    return.params=FALSE,
                                    params = list(),
                                    n.Threads = 4L) {

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

  params <- modifyList(defaults, params)

  ##  Coerce types
  if (!is.matrix(G)) G <- as.matrix(G)
  crosses <- as.matrix(crosses)
  u       <- as.numeric(u)

  if (is.null(fixed.crosses) || nrow(fixed.crosses) == 0) {
    fixed.crosses <- crosses[0, , drop = FALSE]  # 0-row, 2-col
    u.fixed <- u[0]                              # numeric(0)
  }
  fixed.crosses <- as.matrix(fixed.crosses)
  u.fixed <- as.numeric(u.fixed)

  ##  Basic matrix checks
  if (!is.numeric(G)) stop("`G` must be a numeric matrix.")
  if (nrow(G) != ncol(G)) stop("`G` must be square.")
  nInd <- ncol(G)
  if (nInd < 2) stop("`G` must have at least 2 individuals.")

  # Cross matrices must be 2 columns (allow 0 rows for fixed)
  if (ncol(crosses) != 2) stop("`crosses` must have exactly 2 columns.")
  if (nrow(fixed.crosses) > 0L && ncol(fixed.crosses) != 2) {
    stop("`fixed.crosses` must have exactly 2 columns (or 0 rows).")
  }

  ##  Index validity & structure
  if (any(!is.finite(crosses))) stop("`crosses` contains non-finite entries.")
  if (nrow(fixed.crosses) && any(!is.finite(fixed.crosses))) {
    stop("`fixed.crosses` contains non-finite entries.")
  }
  if (any(crosses < 1 | crosses > nInd)) {
    stop("Some indices in `crosses` are outside 1..nrow(G).")
  }
  if (nrow(fixed.crosses) && any(fixed.crosses < 1 | fixed.crosses > nInd)) {
    stop("Some indices in `fixed.crosses` are outside 1..nrow(G).")
  }

  # No self-crosses
  if (any(crosses[,1] == crosses[,2])) stop("`crosses` contains self-crosses.")
  if (nrow(fixed.crosses) && any(fixed.crosses[,1] == fixed.crosses[,2])) {
    stop("`fixed.crosses` contains self-crosses.")
  }

  # Overlap warning
  if (nrow(fixed.crosses)) {
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
  if (nrow(fixed.crosses) != length(u.fixed)) {
    stop("`u.fixed` length (", length(u.fixed), ") must equal nrow(fixed.crosses) (", nrow(fixed.crosses), ").")
  }

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
  if (!is.numeric(params$propability.mutate) || params$propability.mutate < 0 || params$propability.mutate > 1)
    stop("`propability.mutate` must be in [0,1].")
  if (!is.numeric(params$angle.penalty) || params$angle.penalty < 0)
    stop("`angle.penalty` must be >= 0.")

  ##  Cross count logic
  fixedCount <- nrow(fixed.crosses)
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
    Crosses       = crosses,
    fixedCrosses  = fixed.crosses,
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
    res$crossPlan <- as.data.frame(res$crossPlan)[,1:2]
    names(res$crossPlan) <-c("parent1","parent2")
    }else{
    res <- as.data.frame(res$crossPlan)[,1:2]
    }

  return(res)
}
