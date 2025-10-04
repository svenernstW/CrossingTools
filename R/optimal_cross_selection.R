#' Optimal cross selection via a genetic algorithm
#'
#' Runs a GA to select a set of \code{ncrosses} crosses that balance
#' expected gain (\code{g}) and genetic similarity/diversity (via \code{G})
#' around a target trade-off angle.
#'
#'
#' @param ncrosses Integer. Number of crosses to select.
#' @param target.angle Numeric (radians). Target angle between gain and diversity.
#' @param g Numeric vector of per-potential-cross scores (length must be \eqn{nInd*(nInd-1)/2}).
#'   Each element corresponds to one half-diallel cross.
#' @param G Matrix square relationship matrix among individuals (\eqn{nInd \times nInd}).
#' @param propability.mutate Numeric in \code{[0,1]}. Per-solution mutation probability in GA.
#' @param nmutate Integer \(\ge 0\). Number of potential mutations when a mutation occurs.
#' @param nselect Integer \(\ge 2\). Number of solutions retained per generation in GA.
#' @param npop Integer \(\ge 2\). Solution population size (solutions per generation).
#' @param max.generation Integer \(\ge 1\). Maximum GA iterations.
#' @param max.iteration Integer \(\ge 1\). Early-stop if no improvement for this many solution generations.
#' @param angle.penalty Numeric \(\ge 0\). Penalty weight for deviating from \code{target.angle}.
#' @param nthreads Integer \(\ge 1\). number of threads to use.
#'
#' @return A list:
#'   \describe{
#'     \item{crossPlan}{\code{ncrosses x 2} matrix of selected crosses.}
#'     \item{uMax, uMin}{Best/achieved mean cross value}
#'     \item{simMax, simMin}{Best/achieved similarity/relationship.}
#'     \item{uBest, simBest}{Mean cross value and similarity of the final best solution.}
#'     \item{angleBest, lenBest}{Angle and vector length of the final best solution.}
#'   }
#' @examples
#' \dontrun{
#' nInd <- 50
#' G <- diag(nInd) # toy
#' potCross <- nInd*(nInd-1)/2
#' g <- rnorm(potCross)
#' res <- optimal_cross_selection(
#'   ncrosses = 20, target.angle = pi/4, g = g, G = G,
#'   propability.mutate = 0.01, nmutate = 2, nselect = 200,
#'   npop = 2000, max.generation = 200, max.iteration = 50,
#'   angle.penalty = 0.5, nthreads = 4L
#' )
#' }
#' @export
optimal_cross_selection <- function(ncrosses,
                                    target.angle,
                                    g,
                                    G,
                                    propability.mutate = 0.01,
                                    nmutate = 2,
                                    nselect = 500,
                                    npop = 10000,
                                    max.generation = 100,
                                    max.iteration = 100,
                                    angle.penalty,
                                    nthreads = 4L) {
  #  Normalize types (keep your arg names; map to C++ later)
  if (!is.matrix(G)) G <- as.matrix(G)
  g <- as.numeric(g)

  #  Basic shape checks
  # G: square & numeric
  if (!is.numeric(G)) stop("`G` must be a numeric matrix.")
  if (nrow(G) != ncol(G)) stop("`G` must be square: nrow(G) must equal ncol(G).")
  nInd <- ncol(G)
  if (nInd < 2) stop("`G` must have at least 2 individuals (ncol(G) >= 2).")

  # g: length must equal number of potential half-diallel crosses
  potCross <- as.integer(nInd * (nInd - 1) / 2)
  if (length(g) != potCross) {
    stop("`g` must have length nInd*(nInd-1)/2 = ", potCross,
         " (with nInd = ncol(G) = ", nInd, "). Got length ", length(g), ".")
  }
  if (!all(is.finite(g))) stop("`g` contains non-finite values.")

  # Scalars & ranges
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
    if (length(val) != 1L || !is.finite(val)) {
      stop("`", nm, "` must be a single finite value.")
    }
  }

  # Integer checks
  int_params <- c("ncrosses","nmutate","nselect","npop","max.generation","max.iteration","nthreads")
  for (nm in int_params) {
    v <- req_scalar_num[[nm]]
    if (abs(v - round(v)) > .Machine$double.eps^0.5) {
      stop("`", nm, "` must be an integer.")
    }
  }

  ncrosses        <- as.integer(ncrosses)
  nmutate         <- as.integer(nmutate)
  nselect         <- as.integer(nselect)
  npop            <- as.integer(npop)
  max.generation  <- as.integer(max.generation)
  max.iteration   <- as.integer(max.iteration)
  nthreads        <- as.integer(nthreads)

  if (ncrosses < 1L) stop("`ncrosses` must be >= 1.")
  if (ncrosses > potCross) stop("`ncrosses` cannot exceed number of potential crosses (", potCross, ").")
  if (nmutate < 0L) stop("`nmutate` must be >= 0.")
  if (nselect < 2L) stop("`nselect` must be >= 2.")
  if (npop < 2L) stop("`npop` must be >= 2.")
  if (max.generation < 1L) stop("`max.generation` must be >= 1.")
  if (max.iteration < 1L) stop("`max.iteration` must be >= 1.")
  if (nthreads < 1L) stop("`nthreads` must be >= 1.")

  # Probabilities / penalties
  if (!is.numeric(target.angle)) stop("`target.angle` must be numeric.")
  if (!is.numeric(propability.mutate) || propability.mutate < 0 || propability.mutate > 1) {
    stop("`propability.mutate` must be in [0, 1].")
  }
  if (!is.numeric(angle.penalty) || angle.penalty < 0) {
    stop("`angle.penalty` must be >= 0.")
  }


  #Map to C++ parameter names
  nCross       <- ncrosses
  targetAngle  <- as.numeric(target.angle)
  u            <- g
  probMut      <- as.numeric(propability.mutate)
  nMutate      <- nmutate
  nSel         <- nselect
  nPop         <- npop
  maxGen       <- max.generation
  maxRun       <- max.iteration
  anglePenalty <- as.numeric(angle.penalty)
  nThreads     <- nthreads

  #  Call the C++ routine
  res <- cpp_optimal_cross_selection(
    nCross      = nCross,
    targetAngle = targetAngle,
    u           = u,
    G           = G,
    probMut     = probMut,
    nMutate     = nMutate,
    nSel        = nSel,
    nPop        = nPop,
    maxGen      = maxGen,
    maxRun      = maxRun,
    anglePenalty= anglePenalty,
    nThreads    = nThreads
  )


  res
}
