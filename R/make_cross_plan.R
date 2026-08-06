#' Create a crossing plan from candidate parents
#'
#' Generates all possible crosses among a set of candidate parents.
#'
#' For an unsexed design, all unique pairwise crosses are returned, optionally
#' including self-crosses. For a sexed design, all combinations between the
#' supplied male and female parents are returned.
#'
#' @param parents Optional vector of parent identifiers for an unsexed crossing
#'   design. At least two unique identifiers must be supplied.
#' @param self Logical. If \code{TRUE}, include self-crosses when \code{parents}
#'   is supplied. Ignored for sexed crossing designs.
#' @param male.parents Optional vector of male parent identifiers for a sexed
#'   crossing design.
#' @param female.parents Optional vector of female parent identifiers for a
#'   sexed crossing design.
#'
#' @return A data frame with one row per cross. For unsexed designs, the columns
#'   are \code{parent1} and \code{parent2}. For sexed designs, the columns are
#'   \code{male} and \code{female}.
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
