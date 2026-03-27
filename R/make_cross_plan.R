#' Creation of a plan of all potential crosses
#'
#' Generates a crossing plan containing all possible pairwise crosses among a set
#' of candidate genotypes. The function supports both unsexed (symmetric) crossing
#' designs and sexed designs with distinct male and female candidate sets.
#'
#' If \code{parents} is supplied, all unique unordered pairs of parents are returned.
#' Optionally, self-crosses can be included by setting \code{self = TRUE}.
#'
#' Alternatively, if \code{male.parents} and \code{female.parents} are supplied, all possible
#' male–female combinations are returned. In this case, the \code{self} argument
#' is ignored.
#'
#' @param parents A vector of genotype identifiers defining the candidate parents
#'   for an unsexed crossing design. Identifiers may be integers or character strings
#'   (for example, row names of a genotype matrix). Must contain at least two unique
#'   entries.
#'
#' @param self Logical. If \code{TRUE} and \code{parents} is supplied, self-crosses
#'   (parent × itself) are included in addition to all pairwise crosses. Ignored when
#'   \code{male.parents} and \code{female.parents} are supplied.
#'
#' @param male.parents A vector of genotype identifiers defining male parents in a sexed
#'   crossing design. Identifiers may be integers or character strings. Ignored if
#'   \code{parents} is supplied.
#'
#' @param female.parents A vector of genotype identifiers defining female parents in a sexed
#'   crossing design. Identifiers may be integers or character strings. Ignored if
#'   \code{parents} is supplied.
#'
#' @return A two-column \code{data.frame} defining the crossing plan.
#' \itemize{
#'   \item If \code{parents} is supplied, the columns are named \code{parent1} and
#'   \code{parent2} and contain all unique pairwise combinations (and optionally
#'   self-crosses).
#'   \item If \code{male.parents} and \code{female.parents} are supplied, the columns are named
#'   \code{male} and \code{female} and contain all male–female combinations.
#' }
#' The output preserves the type of the supplied identifiers (integer or character),
#' with factors automatically coerced to character.
#'
#' @export

make_cross_plan <- function(
    parents = NULL, self = FALSE, male.parents = NULL, female.parents = NULL
) {
  candidates <- parents
  male_candidates <- male.parents
  female_candidates <- female.parents
  # helper to validate ID vectors
  .validate_ids <- function(x, name, min_len = 1) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.atomic(x)) stop(sprintf("`%s` must be a vector.", name), call. = FALSE)
    if (anyNA(x)) stop(sprintf("`%s` contains NA.", name), call. = FALSE)
    if (length(x) < min_len) {
      stop(sprintf("`%s` must have length >= %d.", name, min_len), call. = FALSE)
    }
    if (any(duplicated(x))) stop(sprintf("`%s` contains duplicated IDs.", name), call. = FALSE)
    invisible(NULL)
  }

  have_candidates <- !is.null(candidates)
  have_male       <- !is.null(male_candidates)
  have_female     <- !is.null(female_candidates)

  if (!have_candidates && !have_male && !have_female) {
    stop("No parents supplied. Provide either `parents` or both `male.parents` and `female.parents`.", call. = FALSE)
  }


  if (have_candidates) {
    if (have_male || have_female) {
      warning("`parents` supplied; ignoring `male.parents` and `female.parents`.")
    }

    if (is.factor(candidates)) candidates <- as.character(candidates)

    .validate_ids(candidates, "parents", min_len = 2)
    ids <- unique(candidates)
    if (length(ids) < 2) stop("`parents` must contain at least 2 unique IDs.", call. = FALSE)
    comb <- t(utils::combn(ids, 2))
    if(self){
      comb <- rbind(comb,cbind(ids,ids))
    }
    plan <- data.frame(
      parent1 = comb[,1],
      parent2 = comb[,2],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    return(plan)
  }


  if (have_male && !have_female) {
    stop("`male.parents` supplied but `female.parents` is missing.", call. = FALSE)
  }
  if (!have_male && have_female) {
    stop("`female.parents` supplied but `male.parents` is missing.", call. = FALSE)
  }
  if(self){
    warning("`female.parents` and `male.parents` supplied, ignoring self argument")
  }

  if (is.factor(male_candidates)) male_candidates <- as.character(male_candidates)
  if (is.factor(female_candidates)) female_candidates <- as.character(female_candidates)

  .validate_ids(male_candidates, "male.parents", min_len = 1)
  .validate_ids(female_candidates, "female.parents", min_len = 1)

  male_ids   <- unique(male_candidates)
  female_ids <- unique(female_candidates)

  plan <- expand.grid(
    male   = male_ids,
    female = female_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(plan) <- NULL
  plan
}
