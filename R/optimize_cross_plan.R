#' Optimal cross selection via a genetic algorithm (with optional fixed crosses)
#'
#' Selects a mating plan consisting of \code{ncrosses} parent pairs by balancing
#' a cross criterion (\code{criterion}) against a genomic similarity (or
#' inbreeding proxy) taken from \code{G.mat}. The optimization can be run in
#' one of two modes:
#' \describe{
#'   \item{\code{method = "angle"}}{Returns a single mating plan that targets a
#'   user-specified trade-off between cross criterion and similarity, expressed as
#'   \code{target.angle} (0--90 degrees). Smaller angles prioritize cross criterion,
#'   larger angles prioritize lower similarity (more diversity).}
#'   \item{\code{method = "pareto"}}{Returns an approximation of the Pareto frontier:
#'   a set of nondominated mating plans where no plan has both (i) lower cross criterion
#'   and (ii) higher similarity than another. This mode is more computationally
#'   demanding but provides multiple optimal trade-offs.}
#' }
#'
#' If \code{fixed.crosses} is supplied, these crosses are always included in every
#' returned plan. The GA then selects the remaining \code{ncrosses - nrow(fixed.crosses)}
#' crosses from \code{crosses}.
#'
#' Both \code{crosses} and \code{fixed.crosses} can be provided either as
#' 1-based integer indices (referring to rows/columns of \code{G.mat}) or as
#' character identifiers matching \code{rownames(G.mat)}.
#'
#' @param candidate.crosses Integer or character matrix/data.frame with two columns (\eqn{K \times 2}).
#'   Candidate (variable) crosses. Must have exactly two columns; self-crosses are not allowed.
#' @param fixed.crosses Integer or character matrix/data.frame with two columns (\eqn{F \times 2}) or \code{NULL}.
#'   Crosses that are always included in the final plan. Can have 0 rows. Must use the
#'   same indexing/identifiers as \code{candidate.crosses}. Self-crosses are not allowed.
#' @param ncrosses Integer. Total number of crosses in the final plan (fixed + variable).
#'   Must be \eqn{\ge} \code{nrow(fixed.crosses)}.
#' @param target.angle Numeric. Target trade-off angle in degrees (0--90) used only when
#'   \code{method = "angle"}. Use \code{0} to prioritize cross criterion only and \code{90} to prioritize
#'   minimizing similarity (maximizing diversity).
#' @param criterion Numeric vector of length \code{nrow(candidate.crosses)}. cross criterion for each candidate cross.
#' @param criterion.fixed Numeric vector of length \code{nrow(fixed.crosses)} (or \code{NULL}).
#'   cross criterion for each fixed cross. If \code{fixed.crosses} has 0 rows, this can be length 0.
#' @param G.mat Numeric square matrix (\eqn{n \times n}). Genomic relationship matrix used to compute
#'   similarity/inbreeding of parental contributions. If \code{candidate.crosses} or \code{fixed.crosses} are
#'   character, \code{rownames(G.mat)} must be present.
#' @param method Character string, either \code{"angle"} or \code{"pareto"}.
#' @param plot Logical. If \code{TRUE} and \code{method = "pareto"}, plot the estimated Pareto frontier.
#'   Ignored when \code{method = "angle"}.
#' @param params List of GA parameters:
#'   \describe{
#'     \item{\code{propability.mutate}}{Numeric in \eqn{[0,1]}. Mutation probability per offspring.}
#'     \item{\code{n.mutate}}{Integer \eqn{\ge 0}. Number of positions to mutate when mutation occurs.}
#'     \item{\code{n.select}}{Integer \eqn{\ge 2}. Number of parents selected each generation.}
#'     \item{\code{n.pop}}{Integer \eqn{\ge 2}. Population size.}
#'     \item{\code{max.generation}}{Integer \eqn{\ge 1}. Maximum number of generations.}
#'     \item{\code{max.iteration}}{Integer \eqn{\ge 1}. Stop after this many generations without improvement.}
#'     \item{\code{angle.penalty}}{Numeric \eqn{\ge 0}. Penalty strength for deviating from \code{target.angle}
#'       (used only when \code{method = "angle"}).}
#'   }
#' @param return.params Logical. If \code{TRUE}, return additional diagnostics (see below).
#' @param nthreads Integer \eqn{\ge 1}. Number of threads to use.
#'
#' @return
#' If \code{return.params = FALSE}:
#' \itemize{
#'   \item A data.frame with columns \code{parent1}, \code{parent2} describing the selected mating plan.
#' }
#'
#' If \code{return.params = TRUE} and \code{method = "angle"}:
#' \itemize{
#'   \item \code{crossPlan}: selected plan (data.frame with \code{parent1}, \code{parent2})
#'   \item \code{uMax}, \code{simMax}: cross criterion optimum and its similarity
#'   \item \code{uBest}, \code{simBest}: achieved solution
#'   \item \code{angleBest}, \code{lenBest}: objective-space angle/length diagnostics
#' }
#'
#' If \code{return.params = TRUE} and \code{method = "pareto"}:
#' \itemize{
#'   \item \code{pareto.plans}: list of mating plans (each a data.frame with \code{parent1}, \code{parent2})
#'   \item \code{pareto.frontier}: data.frame with columns \code{pareto.id}, \code{u}, \code{sim}
#' }
#'
#' @export



optimize_cross_plan <- function(candidate.crosses,
                                    fixed.crosses = NULL,
                                    ncrosses,
                                    target.angle=0,
                                    criterion,
                                    criterion.fixed = NULL,
                                    G.mat,
                                    method = "pareto",
                                    return.params=FALSE,
                                    plot=TRUE,
                                    params = list(),
                                    nthreads = 4L) {

  defaults <-
    list(
      propability.mutate = 0.01,
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
  params <- modifyList(defaults, params)
  u <- criterion
  u.fixed = criterion.fixed
  crosses <- candidate.crosses
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
  if(method=="pareto"){

    ##  Map names to C++
    res <- cpp_optimal_cross_pareto(
      Crosses       = crosses2,
      fixedCrosses  = fixed.crosses2,
      nCross        = ncrosses,
      u             = u,
      ufixed        = u.fixed,
      G             = G,
      probMut       = as.numeric(params$propability.mutate),
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


      paretoPlans  <- lapply (res$paretoPlans, function(x){tempplan <- as.data.frame(x)
                                                            idx <- match(
                                                              paste(tempplan[,1], tempplan[,2]),
                                                              paste(temp[,1],     temp[,2])
                                                            )
                                                            tempplan <- as.data.frame(temp2[idx,])
                                                            names(tempplan) <-c("parent1","parent2")
                                                            return(tempplan)})

      res2 <- list(pareto.plans=paretoPlans, pareto.frontier= pareto)


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


      paretoPlans  <- lapply (res$paretoPlans, function(x){tempplan <- as.data.frame(x)
      idx <- match(
        paste(tempplan[,1], tempplan[,2]),
        paste(temp[,1],     temp[,2])
      )
      tempplan <- as.data.frame(temp2[idx,])
      names(tempplan) <-c("parent1","parent2")
      return(tempplan)})

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
