#' Optimal cross selection via a genetic algorithm (with fixed crosses)
#'
#' @param crosses        integer matrix (K x 2): candidate variable crosses (1-based indices)
#' @param fixed.crosses  integer matrix (F x 2): crosses that must be in the solution (can be 0 rows)
#' @param ncrosses       integer: total crosses in final plan (variable + fixed)
#' @param target.angle   numeric (radians): target tradeoff angle
#' @param u              numeric length K: utility for each row of `crosses`
#' @param u.fixed        numeric length F: utility for each row of `fixed.crosses` (can be length 0)
#' @param G              numeric square GRM (nInd x nInd)
#' @param propability.mutate numeric in \[0-1]
#' @param nmutate, nselect, npop, max.generation, max.iteration integers
#' @param angle.penalty  numeric >= 0
#' @param nthreads       integer >= 1
#' @return list from cpp_optimal_cross_selection()
#' @export
optimal_cross_selection <- function(crosses,
                                    fixed.crosses=NULL,
                                    ncrosses,
                                    target.angle,
                                    u,
                                    u.fixed=NULL,
                                    G,
                                    propability.mutate = 0.01,
                                    nmutate = 2,
                                    nselect = 500,
                                    npop = 10000,
                                    max.generation = 100,
                                    max.iteration = 100,
                                    angle.penalty = 0.5,
                                    nthreads = 4L) {

  ##  Coerce types
  if (!is.matrix(G)) G <- as.matrix(G)
  crosses       <- as.matrix(crosses)
  u       <- as.numeric(u)

  if (is.null(fixed.crosses) || nrow(fixed.crosses) == 0) {
    # keep same column structure (2 cols) and class
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
  # Finite & in-range indices
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



  # Overlap between fixed and variable candidates (warn or stop)
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

  ##  Scalar checks
  req_scalar_num <- list(
    ncrosses = ncrosses,
    target.angle = target.angle,
    propability.mutate = propability.mutate,
    nmutate = nmutate,
    nselect = nselect,
    npop = npop,
    max.generation = max.generation,
    max.iteration = max.iteration,
    angle.penalty = angle.penalty,
    nthreads = nthreads
  )
  for (nm in names(req_scalar_num)) {
    val <- req_scalar_num[[nm]]
    if (length(val) != 1L || !is.finite(val)) stop("`", nm, "` must be a single finite value.")
  }
  # integer-ish
  for (nm in c("ncrosses","nmutate","nselect","npop","max.generation","max.iteration","nthreads")) {
    v <- req_scalar_num[[nm]]
    if (abs(v - round(v)) > .Machine$double.eps^0.5) stop("`", nm, "` must be an integer.")
  }

  ncrosses       <- as.integer(ncrosses)
  nmutate        <- as.integer(nmutate)
  nselect        <- as.integer(nselect)
  npop           <- as.integer(npop)
  max.generation <- as.integer(max.generation)
  max.iteration  <- as.integer(max.iteration)
  nthreads       <- as.integer(nthreads)

  if (ncrosses < 1L) stop("`ncrosses` must be >= 1.")
  if (nmutate < 0L) stop("`nmutate` must be >= 0.")
  if (nselect < 2L) stop("`nselect` must be >= 2.")
  if (npop < 2L) stop("`npop` must be >= 2.")
  if (max.generation < 1L) stop("`max.generation` must be >= 1.")
  if (max.iteration  < 1L) stop("`max.iteration` must be >= 1.")
  if (nthreads < 1L) stop("`nthreads` must be >= 1.")
  if (!is.numeric(target.angle)) stop("`target.angle` must be numeric.")
  if (!is.numeric(propability.mutate) || propability.mutate < 0 || propability.mutate > 1)
    stop("`propability.mutate` must be in [0,1].")
  if (!is.numeric(angle.penalty) || angle.penalty < 0)
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
    probMut       = as.numeric(propability.mutate),
    nMutate       = nmutate,
    nSel          = nselect,
    nPop          = npop,
    maxGen        = max.generation,
    maxRun        = max.iteration,
    anglePenalty  = as.numeric(angle.penalty),
    nThreads      = nthreads
  )

  res
}
