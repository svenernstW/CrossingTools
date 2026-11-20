#' Optimise crossing plans for genetic gain and diversity
#'
#' Selects a crossing plan of \code{ncrosses} crosses by balancing a
#' user-defined cross criterion against genomic similarity among the
#' contributing parents. Fixed crosses and upper and lower constraints on
#' parental contributions or other constraints can be incorporated directly
#' into the optimisation.
#'
#' Ranking crosses solely according to their predicted value may maximise
#' short-term genetic gain but can increase coancestry among selected parents
#' and reduce genetic diversity. Optimal cross selection therefore combines
#' expected cross performance with a measure of relatedness among the parents
#' contributing to a crossing plan.
#'
#' For each candidate crossing plan, performance is quantified by the mean
#' value of the supplied cross criterion across all selected crosses. Genomic
#' similarity is calculated from proportional parental contributions and the
#' relationship matrix supplied in \code{G.mat}, analogous to group coancestry
#' as used in optimal contribution selection (Meuwissen, 1997). Each mating
#' contributes equally to the crossing plan, with half of its contribution
#' assigned to each parent.
#'
#' Because the cross criterion and genomic similarity may be measured on
#' different scales, the objectives are normalised for the angle-based
#' optimisation. CrossingTools provides two strategies for balancing the
#' gain--diversity objectives. With \code{method = "angle"}, the user specifies
#' a target angle between 0 and 90 degrees. An angle of 0 degrees gives maximum
#' emphasis to the cross criterion, whereas an angle of 90 degrees gives
#' maximum emphasis to reducing genomic similarity and therefore maintaining
#' diversity. Intermediate angles specify different compromises between the
#' two objectives.
#'
#' With \code{method = "pareto"}, the cross criterion and genomic similarity
#' are treated as separate objectives and optimised simultaneously using an
#' NSGA-II-type multi-objective genetic algorithm (Deb et al., 2002). The
#' returned non-dominated solutions approximate the Pareto frontier and
#' represent alternative trade-offs between predicted cross performance and
#' genomic diversity.
#'
#' Both optimisation strategies search over complete crossing plans rather
#' than ranking crosses independently. Practical breeding constraints can be
#' imposed through fixed crosses and minimum or maximum parental contributions.
#' Additional constraints can be specified through \code{cross.constraints},
#' allowing the mean of numerical cross properties or the proportion of
#' individual cross types to be restricted to user-defined intervals.
#' These constraints can be applied simultaneously with parental contribution
#' limits and are treated as hard constraints rather than additional
#' optimisation objectives.
#'
#' Numerical constraints can represent any quantitative cross property, such
#' as expected trait values, expected allele frequencies of QTL, genetic
#' distances, or crossing costs. Categorical constraints can be used to restrict
#' the proportions of particular cross types, such as within-family or
#' between-family crosses. Individual crosses may have numerical values outside
#' the specified bounds as long as the overall mean satisfies the constraint.
#'
#' Crosses that should not be considered can be excluded from
#' \code{candidate.crosses}.
#'
#' @param candidate.crosses Integer or character matrix or data frame with two
#'   columns. Candidate crosses from which the variable part of the crossing
#'   plan is selected. Self-crosses are not allowed.
#' @param fixed.crosses Integer or character matrix or data frame with two
#'   columns, or \code{NULL}. Crosses included in every returned plan. These
#'   must use the same type of parent identifiers as \code{candidate.crosses}.
#'   Self-crosses are not allowed.
#' @param ncrosses Integer. Total number of crosses in each returned plan,
#'   including fixed crosses. Must be at least \code{nrow(fixed.crosses)}.
#' @param target.angle Numeric scalar between 0 and 90. Target trade-off angle
#'   used when \code{method = "angle"}. An angle of 0 gives maximum emphasis
#'   to the cross criterion, whereas an angle of 90 gives maximum emphasis to
#'   reducing genomic similarity.
#' @param criterion Numeric vector with one cross-criterion value for each row
#'   of \code{candidate.crosses}. Larger values are considered favourable.
#' @param criterion.fixed Numeric vector with one cross-criterion value for each
#'   row of \code{fixed.crosses}, or \code{NULL} when no fixed crosses are used.
#' @param G.mat Numeric square matrix containing genomic relationships,
#'   similarities, or another pairwise measure used to quantify similarity
#'   among parental contributions. Character parent identifiers require
#'   matching row names.
#' @param parents.upper Numeric vector with one value per parent in
#'   \code{G.mat}, giving the maximum number of times each parent may occur
#'   across the complete crossing plan. Each cross contributes one occurrence
#'   to each parent. Values must be non-negative whole numbers or \code{Inf}.
#'   The default is \code{Inf} for every parent.
#' @param parents.lower Numeric vector with one value per parent in
#'   \code{G.mat}, giving the minimum number of times each parent must occur
#'   across the complete crossing plan. Each cross contributes one occurrence
#'   to each parent. Values must be non-negative whole numbers and cannot exceed
#'   the corresponding values in \code{parents.upper}. The default is zero for
#'   every parent.
#' @param method Character string specifying the optimisation method. Must be
#'   either \code{"angle"} or \code{"pareto"}.
#' @param cross.constraints Optional named list specifying additional
#'   constraints on numerical or categorical properties of the selected
#'   crosses. The default is \code{NULL}, corresponding to no additional
#'   cross-level constraints. The list may contain the following elements:
#'   \describe{
#'     \item{\code{values}}{
#'       A data frame or matrix with one row per candidate cross and one
#'       column per cross property. Rows must correspond to the order of
#'       \code{candidate.crosses}. Columns may contain numeric values,
#'       character strings, or factors. Numerical columns are interpreted
#'       as quantitative cross properties, whereas character or factor
#'       columns define categorical cross types.
#'     }
#'     \item{\code{fixed.values}}{
#'       A data frame or matrix containing the corresponding properties
#'       of fixed crosses, with one row per cross in \code{fixed.crosses}.
#'       The columns must match those in \code{values}. This element is
#'       required when fixed crosses are supplied and may otherwise be
#'       omitted or set to \code{NULL}.
#'     }
#'     \item{\code{bounds}}{
#'       A named list specifying the constraints to be imposed.
#'       Each list name must correspond to a column in \code{values}.
#'       The required format depends on whether the corresponding
#'       property is numerical or categorical.
#'     }
#'   }
#'
#'   For numerical properties, each element of \code{bounds} must be a
#'   named numeric vector of the form
#'   \code{c(lower = ..., upper = ...)}, specifying the permissible
#'   interval for the mean property value across the complete crossing
#'   plan. Bounds are expressed in the original units of the property.
#'   The values \code{-Inf} and \code{Inf} can be used to specify
#'   unrestricted lower or upper bounds, respectively.
#'
#'   For categorical properties, each element of \code{bounds} must be a
#'   data frame with columns \code{type}, \code{lower}, and \code{upper}.
#'   Each row specifies a category and its permissible proportion among
#'   selected crosses. Category bounds must lie between zero and one.
#'   Categories not included in the bounds table are not constrained
#'   directly.
#'
#'   Numerical and categorical constraints may be combined.
#'   Only properties explicitly included in \code{bounds} are constrained.
#'   Constraints are applied to the complete crossing plan, including
#'   fixed crosses. The user is responsible for ensuring that each row
#'   of \code{values} and \code{fixed.values} corresponds to the correct
#'   cross. All additional constraints must be jointly feasible for
#'   a valid crossing plan to exist.
#' @param return.params Logical. Whether to return optimisation diagnostics in
#'   addition to the selected crossing plan or Pareto solutions.
#' @param plot Logical. Whether to display the estimated Pareto frontier when
#'   \code{method = "pareto"}. Ignored when \code{method = "angle"}.
#' @param params Named list of genetic-algorithm parameters. Unspecified entries
#'   use their default values:
#'   \describe{
#'     \item{\code{probability.mutate}}{Numeric scalar in \eqn{[0,1]}. Mutation
#'     probability for each offspring.}
#'     \item{\code{n.mutate}}{Non-negative integer. Number of positions changed
#'     when mutation occurs.}
#'     \item{\code{n.select}}{Integer of at least 2. Number of solutions selected
#'     in each generation.}
#'     \item{\code{n.pop}}{Integer of at least 2. Population size.}
#'     \item{\code{max.generation}}{Positive integer. Maximum number of
#'     generations.}
#'     \item{\code{max.iteration}}{Positive integer. Number of generations
#'     without improvement after which optimisation stops.}
#'     \item{\code{angle.penalty}}{Non-negative numeric scalar controlling the
#'     penalty for deviation from \code{target.angle}, used only when
#'     \code{method = "angle"}.}
#'   }
#' @param nthreads Positive integer. Number of computational threads.
#'
#' @return
#' For \code{method = "angle"} and \code{return.params = TRUE}, a list
#' containing:
#' \describe{
#'   \item{\code{crossPlan}}{The selected crossing plan as a data frame with
#'   columns \code{parent1} and \code{parent2}.}
#'   \item{\code{uMax}, \code{simMax}}{The cross criterion of the
#'   high-gain reference plan and its corresponding genomic similarity.
#'   With active constraints, the reference plan is constructed
#'   heuristically and is not guaranteed to be the global gain maximum.}
#'   \item{\code{uMin}, \code{simMin}}{Reference values obtained during
#'   similarity optimisation, used for normalisation of the objectives.}
#'   \item{\code{uBest}, \code{simBest}}{The cross criterion and similarity of
#'   the selected solution.}
#'   \item{\code{angleBest}, \code{lenBest}}{The angle and length of the
#'   selected solution in normalised objective space.}
#'   \item{\code{parentContributions}}{Numeric vector containing the
#'   proportional contribution of each parent to the selected crossing
#'   plan, including fixed crosses. Contributions sum to one.
#'   Parent names are retained when available in \code{G.mat}.}
#'   \item{\code{constraintValues}}{Named numeric vector containing the
#'   achieved mean values of numerical constraints and the proportions
#'   of constrained categorical cross types. Numerical values are
#'   reported in their original units. An empty vector is returned when
#'   no additional constraints are specified.}
#' }
#'
#' For \code{target.angle = 0}, the high-gain reference plan is returned
#' directly without running the genetic algorithm. In this case,
#' \code{uMax}, \code{simMax}, \code{parentContributions}, and
#' \code{constraintValues} are returned with \code{crossPlan};
#' the remaining optimisation diagnostics are not available.
#'
#' For \code{method = "pareto"}, a list containing:
#' \describe{
#'   \item{\code{pareto.plans}}{A list of non-dominated crossing plans
#'   identified by the optimisation, each represented by a data frame
#'   with columns \code{parent1} and \code{parent2}. All returned plans
#'   satisfy the specified parental and additional cross constraints.}
#'   \item{\code{pareto.frontier}}{A data frame containing the solution
#'   identifier \code{pareto.id}, cross criterion \code{u}, and genomic
#'   similarity \code{sim}. Each row corresponds to the crossing plan
#'   at the same position in \code{pareto.plans}.}
#'   \item{\code{pareto.constraints}}{When \code{return.params = TRUE} and
#'   additional cross constraints are specified, a data frame containing
#'   the achieved constraint values for every returned Pareto solution.
#'   Rows correspond to \code{pareto.frontier} and \code{pareto.plans},
#'   while columns represent individual numerical constraints or
#'   constrained categorical cross types. Numerical properties are
#'   reported in their original units and categorical properties as
#'   proportions. This element is omitted when no additional cross
#'   constraints are specified.}
#' }
#'
#' @references
#' Meuwissen, T. H. E. (1997).
#' Maximizing the response of selection with a predefined rate of inbreeding.
#' \emph{Journal of Animal Science}, 75(4), 934--940.
#'
#' Deb, K., Pratap, A., Agarwal, S. and Meyarivan, T. (2002).
#' A fast and elitist multiobjective genetic algorithm: NSGA-II.
#' \emph{IEEE Transactions on Evolutionary Computation}, 6(2), 182--197.
#' \doi{10.1109/4235.996017}
#'
#' @export






optimize_cross_plan <- function(candidate.crosses,
                                fixed.crosses = NULL,
                                ncrosses,
                                target.angle=0,
                                criterion,
                                criterion.fixed = NULL,
                                G.mat,
                                parents.upper=NA,
                                parents.lower=NA,
                                method = "pareto",
                                return.params=FALSE,
                                plot=TRUE,
                                params = list(),
                                nthreads = 4L,
                                cross.constraints = NULL) {

  defaults <-
    list(
      probability.mutate = 0.01,
      n.mutate = 2,
      n.select = 400,
      n.pop = 2000,
      max.generation = 100,
      max.iteration = 20,
      angle.penalty = 0.5
    )


  if(!method %in% c("pareto","angle")){
    stop("Method needs to be pareto or angle!")
  }

  G <- G.mat
  n.Threads <- nthreads
  params <- utils::modifyList(defaults, params)
  u <- criterion
  u.fixed = criterion.fixed
  crosses <- candidate.crosses
  ##  Coerce types

  nInd <- nrow(G)

  # Default: no lower or upper constraints
  if (length(parents.upper) == 1L && is.na(parents.upper)) {
    parents.upper <- rep(Inf, nInd)
  }

  if (length(parents.lower) == 1L && is.na(parents.lower)) {
    parents.lower <- rep(0, nInd)
  }

  if (!is.numeric(parents.upper) || !is.numeric(parents.lower)) {
    stop("`parents.upper` and `parents.lower` must be numeric.")
  }

  if (length(parents.upper) != nInd) {
    stop(
      "`parents.upper` must have exactly one value per parent in `G.mat`."
    )
  }

  if (length(parents.lower) != nInd) {
    stop(
      "`parents.lower` must have exactly one value per parent in `G.mat`."
    )
  }

  if (anyNA(parents.upper) || anyNA(parents.lower)) {
    stop(
      "`parents.upper` and `parents.lower` cannot contain NA or NaN."
    )
  }

  if (any(!is.finite(parents.lower))) {
    stop("All values in `parents.lower` must be finite.")
  }

  if (any(parents.lower < 0) || any(parents.upper < 0)) {
    stop("Parent contribution bounds cannot be negative.")
  }

  if (any(parents.lower > parents.upper)) {
    stop(
      "Every value in `parents.lower` must be <= the corresponding ",
      "value in `parents.upper`."
    )
  }

  # The bounds represent numbers of parental appearances and must
  # therefore be whole numbers, except that upper bounds may be Inf.
  if (any(abs(parents.lower - round(parents.lower)) >
          .Machine$double.eps^0.5)) {
    stop("All values in `parents.lower` must be whole numbers.")
  }

  finiteUpper <- is.finite(parents.upper)

  if (any(
    abs(parents.upper[finiteUpper] -
        round(parents.upper[finiteUpper])) >
    .Machine$double.eps^0.5
  )) {
    stop(
      "All finite values in `parents.upper` must be whole numbers."
    )
  }

  parents.lower <- as.numeric(round(parents.lower))
  parents.upper[finiteUpper] <-
    round(parents.upper[finiteUpper])

  # There are exactly 2 * ncrosses parental appearances.
  if (sum(parents.lower) > 2 * ncrosses) {
    stop(
      "The specified `parents.lower` values require more parental ",
      "appearances than are available in the requested mating plan."
    )
  }

  if (sum(parents.upper) < 2 * ncrosses) {
    stop(
      "The specified `parents.upper` values permit fewer parental ",
      "appearances than are required in the requested mating plan."
    )
  }







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
  if (nrow(fixed.crosses2) > 0 && length(u.fixed) != nrow(fixed.crosses2)) {
    stop("`criterion.fixed` must have length nrow(fixed.crosses).")
  }

  u.fixed <- as.numeric(u.fixed)

  ##  Basic matrix checks
  if (!is.numeric(G)) stop("`G.mat` must be a numeric matrix.")
  if (nrow(G) != ncol(G)) stop("`G.mat` must be square.")
  nInd <- ncol(G)
  if (nInd < 2) stop("`G.mat` must have at least 2 individuals.")

  # Cross matrices must be 2 columns (allow 0 rows for fixed)
  if (ncol(crosses) != 2) stop("`crosses` must have exactly 2 columns.")
  if (!is.null(fixed.crosses) && nrow(fixed.crosses) > 0L && ncol(fixed.crosses) != 2) {
    stop("`fixed.crosses` must have exactly 2 columns (or 0 rows).")
  }


  ##  Index validity & structure
  if (any(!is.finite(crosses2))) stop("`crosses` contains non-finite entries.")
  if (nrow(fixed.crosses2) && any(!is.finite(fixed.crosses2))) stop("`fixed.crosses` contains non-finite entries.")

  if (any(crosses2 < 1 | crosses2 > nInd)) stop("Some indices in `crosses` are outside 1..nrow(G.mat).")
  if (nrow(fixed.crosses2) && any(fixed.crosses2 < 1 | fixed.crosses2 > nInd)) {
    stop("Some indices in `fixed.crosses` are outside 1..nrow(G.mat).")
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
    stop("`criterion` length (", length(u), ") must equal nrow(crosses) (", nrow(crosses), ").")
  }
  if (ncol(crosses) != 2) stop("`crosses` must have exactly 2 columns.")
  if (nrow(fixed.crosses2) > 0L && ncol(fixed.crosses2) != 2) stop("`fixed.crosses` must have exactly 2 columns (or 0 rows).")


  ##  Scalar checks (use n.Threads)
  req_scalar_num <- list(
    ncrosses = ncrosses,
    target.angle = target.angle,
    probability.mutate = params$probability.mutate,
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
  if (!is.numeric(params$probability.mutate) || params$probability.mutate < 0 || params$probability.mutate > 1){
    stop("`probability.mutate` must be in [0,1].")}
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


  fixedCount <- nrow(fixed.crosses2)
  nVar <- ncrosses - fixedCount
  fixedParentCount <- tabulate(
    as.vector(fixed.crosses2),
    nbins = nInd
  )

  names(fixedParentCount) <- rownames(G)


  minContr <- parents.lower / (2 * ncrosses)
  maxContr <- parents.upper / (2 * ncrosses)

  # Fixed crosses alone must not exceed an upper bound.
  fixedUpperViolation <-
    fixedParentCount > parents.upper

  if (any(fixedUpperViolation)) {
    bad <- which(fixedUpperViolation)

    badNames <- if (!is.null(rownames(G))) {
      rownames(G)[bad]
    } else {
      bad
    }

    stop(
      "The fixed crosses already exceed `parents.upper` for: ",
      paste(badNames, collapse = ", "),
      ". Fixed parental counts: ",
      paste(fixedParentCount[bad], collapse = ", "),
      "; upper bounds: ",
      paste(parents.upper[bad], collapse = ", "),
      "."
    )
  }

  #############################################################################
  # Additional cross constraints
  #############################################################################

  # Initialize empty constraint matrices.
  # These remain empty when cross.constraints = NULL.

  A <- matrix(numeric(0), nrow = nrow(crosses2), ncol = 0)
  Afixed <- matrix(numeric(0), nrow = nrow(fixed.crosses2), ncol = 0)

  lower <- numeric(0)
  upper <- numeric(0)

  constraint.names <- character(0)

  # Store the original scales so that numerical constraint
  # values can be returned in their original units.
  constraint.offset <- numeric(0)
  constraint.scale <- numeric(0)


  if (!is.null(cross.constraints)) {

    ###########################################################################
    # Validate the list structure
    ###########################################################################

    if (
      !is.list(cross.constraints) ||
      !all(c("values", "bounds") %in% names(cross.constraints))
    ) {
      stop(
        "`cross.constraints` must be a list containing ",
        "`values` and `bounds`."
      )
    }

    values <- cross.constraints$values
    fixed.values <- cross.constraints$fixed.values
    bounds <- cross.constraints$bounds

    if (!is.data.frame(values) && !is.matrix(values)) {
      stop("`cross.constraints$values` must be a data frame or matrix.")
    }

    values <- as.data.frame(values, stringsAsFactors = FALSE)

    if (
      nrow(values) != nrow(crosses2) ||
      is.null(names(values)) ||
      anyNA(names(values)) ||
      any(!nzchar(names(values))) ||
      anyDuplicated(names(values))
    ) {
      stop(
        "`cross.constraints$values` must have one row per candidate ",
        "cross and unique, non-empty column names."
      )
    }

    if (
      !is.list(bounds) ||
      is.data.frame(bounds) ||
      (length(bounds) > 0 &&
       (is.null(names(bounds)) ||
        anyNA(names(bounds)) ||
        any(!nzchar(names(bounds))) ||
        anyDuplicated(names(bounds)))) ||
      !all(names(bounds) %in% names(values))
    ) {
      stop(
        "`cross.constraints$bounds` must be a named list ",
        "corresponding to columns in `values`."
      )
    }


    ###########################################################################
    # Handle fixed crosses
    ###########################################################################

    if (is.null(fixed.values)) {

      if (fixedCount > 0L) {
        stop(
          "`cross.constraints$fixed.values` must be supplied ",
          "when fixed crosses are present."
        )
      }

      fixed.values <- values[0, , drop = FALSE]
    }

    if (!is.data.frame(fixed.values) && !is.matrix(fixed.values)) {
      stop(
        "`cross.constraints$fixed.values` must be a data frame or matrix."
      )
    }

    fixed.values <- as.data.frame(
      fixed.values,
      stringsAsFactors = FALSE
    )

    if (
      nrow(fixed.values) != fixedCount ||
      !setequal(names(fixed.values), names(values)) ||
      anyDuplicated(names(fixed.values))
    ) {
      stop(
        "`cross.constraints$fixed.values` must have one row per ",
        "fixed cross and the same columns as `values`."
      )
    }

    # Match the order of the candidate constraint columns.
    fixed.values <- fixed.values[
      ,
      names(values),
      drop = FALSE
    ]


    ###########################################################################
    # Process every bounded property
    ###########################################################################

    for (nm in names(bounds)) {

      x <- values[[nm]]
      xf <- fixed.values[[nm]]
      b <- bounds[[nm]]


      #########################################################################
      # Numeric constraints
      #########################################################################

      if (is.numeric(x) && !is.factor(x)) {

        if (
          !is.numeric(xf) ||
          anyNA(x) ||
          anyNA(xf) ||
          any(!is.finite(x)) ||
          any(!is.finite(xf))
        ) {
          stop(
            "Numeric constraint `", nm,
            "` must contain finite numeric values."
          )
        }

        # Numeric constraints require a named vector:
        # c(lower = ..., upper = ...)

        if (
          !is.numeric(b) ||
          length(b) != 2L ||
          !setequal(names(b), c("lower", "upper")) ||
          anyNA(b)
        ) {
          stop(
            "Bounds for numeric constraint `", nm,
            "` must be c(lower = ..., upper = ...)."
          )
        }

        lo <- unname(b["lower"])
        hi <- unname(b["upper"])

        if (
          lo > hi ||
          lo == Inf ||
          hi == -Inf
        ) {
          stop(
            "Invalid bounds for numeric constraint `", nm, "`."
          )
        }

        # Use both candidate and fixed crosses to determine
        # the minimum and maximum of the property.

        all.x <- c(x, xf)

        xmin <- min(all.x)
        xmax <- max(all.x)
        span <- xmax - xmin

        if (!is.finite(span)) {
          stop(
            "Numeric range is too large for constraint `", nm, "`."
          )
        }

        # An average of finite values cannot be below xmin or above xmax.
        # These checks are only necessary conditions for feasibility.

        if (lo > xmax || hi < xmin) {
          stop(
            "Numeric constraint `", nm, "` is impossible: ",
            "the available values range from ", xmin,
            " to ", xmax, "."
          )
        }

        if (span == 0) {

          # All crosses have the same value.
          # The bound was already checked above.
          # Use a nonbinding column but retain its name and original value.

          xn <- rep(0, length(x))
          xfn <- rep(0, length(xf))

          ln <- 0
          un <- 1

        } else {

          # Normalize each value to the interval [0,1].

          xn <- (x - xmin) / span
          xfn <- (xf - xmin) / span

          # Apply the same linear transformation to the bounds.
          # Clamp nonbinding bounds to the attainable score range.

          ln <- if (is.infinite(lo)) {
            0
          } else {
            max(0, (lo - xmin) / span)
          }

          un <- if (is.infinite(hi)) {
            1
          } else {
            min(1, (hi - xmin) / span)
          }
        }

        # Append the normalized property to the C++ matrices.

        A <- cbind(
          A,
          matrix(xn, nrow = nrow(crosses2), ncol = 1)
        )

        Afixed <- cbind(
          Afixed,
          matrix(xfn, nrow = fixedCount, ncol = 1)
        )

        lower <- c(lower, ln)
        upper <- c(upper, un)

        constraint.names <- c(constraint.names, nm)

        # Needed to transform the result back to the original scale.
        constraint.offset <- c(constraint.offset, xmin)
        constraint.scale <- c(constraint.scale, span)


        #########################################################################
        # Categorical constraints
        #########################################################################

      } else if (is.character(x) || is.factor(x)) {

        if (
          !(is.character(xf) || is.factor(xf)) ||
          anyNA(x) ||
          anyNA(xf) ||
          any(!nzchar(as.character(x))) ||
          any(!nzchar(as.character(xf)))
        ) {
          stop(
            "Categorical constraint `", nm,
            "` contains missing or invalid categories."
          )
        }

        # Categorical constraints require a data frame
        # containing type, lower and upper.

        if (
          !is.data.frame(b) ||
          !all(c("type", "lower", "upper") %in% names(b))
        ) {
          stop(
            "Bounds for categorical constraint `", nm,
            "` must be a data frame with columns ",
            "`type`, `lower` and `upper`."
          )
        }

        types <- as.character(b$type)

        if (
          length(types) == 0L ||
          anyNA(types) ||
          any(!nzchar(types)) ||
          anyDuplicated(types) ||
          !is.numeric(b$lower) ||
          !is.numeric(b$upper) ||
          anyNA(b$lower) ||
          anyNA(b$upper) ||
          any(!is.finite(b$lower)) ||
          any(!is.finite(b$upper)) ||
          any(b$lower < 0 | b$lower > 1) ||
          any(b$upper < 0 | b$upper > 1) ||
          any(b$lower > b$upper)
        ) {
          stop(
            "Invalid category bounds for constraint `", nm, "`."
          )
        }

        # The categories within a classification are mutually exclusive.
        # Therefore their lower bounds cannot sum to more than 1.

        if (sum(b$lower) > 1 + 1e-12) {
          stop(
            "Lower bounds for categorical constraint `", nm,
            "` sum to more than one."
          )
        }

        # If all observed categories have an upper bound,
        # their permitted proportions must sum to at least 1.

        observed <- unique(c(as.character(x), as.character(xf)))

        if (all(observed %in% types)) {

          if (sum(b$upper[match(observed, types)]) < 1 - 1e-12) {
            stop(
              "Upper bounds for categorical constraint `", nm,
              "` sum to less than one."
            )
          }
        }

        # Convert each constrained category into its own indicator column.

        for (k in seq_along(types)) {

          category <- types[k]

          xn <- as.numeric(as.character(x) == category)
          xfn <- as.numeric(as.character(xf) == category)

          A <- cbind(
            A,
            matrix(xn, nrow = nrow(crosses2), ncol = 1)
          )

          Afixed <- cbind(
            Afixed,
            matrix(xfn, nrow = fixedCount, ncol = 1)
          )

          lower <- c(lower, b$lower[k])
          upper <- c(upper, b$upper[k])

          constraint.names <- c(
            constraint.names,
            paste0(nm, ".", category)
          )

          # Categorical indicators already use the original 0/1 scale.
          constraint.offset <- c(constraint.offset, 0)
          constraint.scale <- c(constraint.scale, 1)
        }

      } else {

        stop(
          "Constraint `", nm,
          "` must be numeric, character or factor."
        )
      }
    }
  }


  #############################################################################
  # Finalize constraint names
  #############################################################################

  constraint.names <- make.unique(constraint.names)

  colnames(A) <- constraint.names
  colnames(Afixed) <- constraint.names

  names(lower) <- constraint.names
  names(upper) <- constraint.names


  if(method=="angle"){

    ##  Map names to C++
    res <- cpp_optimal_cross_selection(
      Crosses       = crosses2,
      fixedCrosses  = fixed.crosses2,
      nCross        = ncrosses,
      targetAngle   = as.numeric(target.angle),
      u             = u,
      ufixed        = u.fixed,
      G             = G,
      minContr      = minContr,
      maxContr      = maxContr,
      A             = A,
      Afixed        = Afixed,
      lower         = lower,
      upper         = upper,
      probMut       = as.numeric(params$probability.mutate),
      nMutate       = nmutate,
      nSel          = nselect,
      nPop          = npop,
      maxGen        = max.generation,
      maxRun        = max.iteration,
      anglePenalty  = as.numeric(params$angle.penalty),
      nThreads      = nThreads
    )

    if(return.params){

      # Give parental contributions their original parent names.
      if (!is.null(rownames(G))) {
        names(res$parentContributions) <- rownames(G)
      }

      # Convert numerical constraints back to their original units.
      # Categorical proportions already have offset 0 and scale 1.
      res$constraintValues <-
        res$constraintValues * constraint.scale + constraint.offset

      names(res$constraintValues) <- constraint.names

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
  if(method=="pareto"){

    ##  Map names to C++
    res <- cpp_optimal_cross_pareto(
      Crosses       = crosses2,
      fixedCrosses  = fixed.crosses2,
      nCross        = ncrosses,
      u             = u,
      ufixed        = u.fixed,
      G             = G,
      minContr      = minContr,
      maxContr      = maxContr,
      A             = A,
      Afixed        = Afixed,
      lower         = lower,
      upper         = upper,
      probMut       = as.numeric(params$probability.mutate),
      nMutate       = nmutate,
      nSel          = nselect,
      nPop          = npop,
      maxGen        = max.generation,
      maxRun        = max.iteration,
      nThreads      = nThreads
    )
    if(return.params){
      paretoPlans  <- lapply (res$paretoPlans, function(x){temp <- as.data.frame(x)
      names(temp)<-c("parent1","parent2")

      return(temp)})

      temp <- rbind(crosses2,fixed.crosses2)
      if (is.null(fixed.crosses) || nrow(fixed.crosses)==0){
        fixed.orig <- crosses[0,,drop=FALSE]
      } else{
        fixed.orig <- fixed.crosses
      }
      temp2 <- rbind(crosses, fixed.orig)

      pareto <- as.data.frame(res$pareto)
      names(pareto)[1] <- "pareto.id"


      labs <- rownames(G.mat)

      paretoPlans <- lapply(res$paretoPlans, function(x) {
        df <- as.data.frame(x)
        names(df) <- c("parent1", "parent2")

        # Only if original input was character
        if (!is.numeric(candidate.crosses) && !is.integer(candidate.crosses)) {
          df$parent1 <- labs[df$parent1]
          df$parent2 <- labs[df$parent2]
        }

        df
      })

      res2 <- list(pareto.plans=paretoPlans, pareto.frontier= pareto)

      # Convert Pareto constraint values back to their original units.

      if (length(constraint.names) > 0L) {

        pareto.constraints <- as.matrix(
          res$paretoConstraintValues
        )

        pareto.constraints <- sweep(
          pareto.constraints,
          2,
          constraint.scale,
          "*"
        )

        pareto.constraints <- sweep(
          pareto.constraints,
          2,
          constraint.offset,
          "+"
        )

        colnames(pareto.constraints) <- constraint.names

        res2$pareto.constraints <- as.data.frame(
          pareto.constraints
        )
      }

    }else{
      paretoPlans  <- lapply (res$paretoPlans, function(x){temp <- as.data.frame(x)
      names(temp)<-c("parent1","parent2")

      return(temp)})

      temp <- rbind(crosses2,fixed.crosses2)
      if (is.null(fixed.crosses) || nrow(fixed.crosses)==0){
        fixed.orig <- crosses[0,,drop=FALSE]
      } else{
        fixed.orig <- fixed.crosses
      }
      temp2 <- rbind(crosses, fixed.orig)

      pareto <- as.data.frame(res$pareto)
      names(pareto)[1] <- "pareto.id"


      key_temp <- paste(temp[,1], temp[,2], sep = "_")

      paretoPlans <- lapply(res$paretoPlans, function(x) {

        x <- as.data.frame(x)
        names(x) <- c("parent1", "parent2")

        key_x <- paste(x$parent1, x$parent2, sep = "_")

        idx <- match(key_x, key_temp)

        if (anyNA(idx)) {
          key_x_rev <- paste(x$parent2, x$parent1, sep = "_")
          idx2 <- match(key_x_rev[is.na(idx)], key_temp)
          idx[is.na(idx)] <- idx2
        }

        out <- as.data.frame(temp2[idx, , drop = FALSE])
        names(out) <- c("parent1", "parent2")
        out
      })
      res2 <- list(pareto.plans=paretoPlans, pareto.frontier= pareto)
    }

    if(plot){
      df <- res2$pareto.frontier
      df <- df[order(df$sim, df$u), ]

      df$label <- paste0(
        "id: ", df$pareto.id,
        "<br>sim: ", signif(df$sim, 5),
        "<br>u: ", signif(df$u, 5)
      )

      p <- ggplot2::ggplot(df, ggplot2::aes(x = sim, y = u)) +
        ggplot2::geom_point(
          ggplot2::aes(text = label),
          size = 1.5, alpha = 0.9, shape = 4
        ) +
        ggplot2::geom_path(ggplot2::aes(group = 1), linewidth = 0.8) +
        ggplot2::labs(
          x = "'inbreeding' (lower = better)",
          y = "u (higher = better)"
        ) +
        ggplot2::theme_grey(base_size = 10) +
        ggplot2::theme(
          legend.position   = "bottom",
          legend.key        = ggplot2::element_rect(fill = "transparent", colour = NA),
          legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
          panel.background  = ggplot2::element_rect(colour = "black", fill = "grey93", linewidth = 1.1),
          axis.title        = ggplot2::element_text(size = 11),
          axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
          axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
          axis.ticks        = ggplot2::element_line(),
          axis.text         = ggplot2::element_text(size = 10)
        )

      print(plotly::ggplotly(p, tooltip = "label"))


    }

    return(res2)


  }

}
