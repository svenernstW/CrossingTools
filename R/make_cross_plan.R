#' Creation of a plan of all potential crosses for CrossingTools
#'
#'
#'
#' @param candidates A vector of integers that identify each genotype. Should be consistent with genotypic data.
#' Use this function if there are no sexes.
#' @param male_candidates A vector of integers that identify each male genotype. Ignored if candidates is supplied.
#' Should be consistent with genotypic data.
#' @param female_candidates A vector of integers that identify each female genotype. Ignored if candidates is supplied.
#' Should be consistent with genotypic data.
#' @return A two-column data.frame with the the combinations of all supplied individuals.
#' If male and female candidates are supplied, the first column corresponds to males, and the second to females
#' @export
make_cross_plan <- function(
    candidates = NULL, male_candidates = NULL, female_candidates = NULL
) {
  # helper to validate ID vectors
  .validate_ids <- function(x, name, min_len = 1) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.atomic(x)) stop(sprintf("`%s` must be an atomic vector.", name), call. = FALSE)
    if (!is.numeric(x) && !is.integer(x)) {
      stop(sprintf("`%s` must be numeric/integer IDs.", name), call. = FALSE)
    }
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
    stop("No candidates supplied. Provide either `candidates` or both `male_candidates` and `female_candidates`.", call. = FALSE)
  }


  if (have_candidates) {
    if (have_male || have_female) {
      warning("`candidates` supplied; ignoring `male_candidates` and `female_candidates`.")
    }
    .validate_ids(candidates, "candidates", min_len = 2)
    ids <- unique(candidates)
    if (length(ids) < 2) stop("`candidates` must contain at least 2 unique IDs.", call. = FALSE)
    comb <- utils::combn(ids, 2)
    plan <- data.frame(
      parent1 = comb[1, ],
      parent2 = comb[2, ],
      row.names = NULL
    )
    return(plan)
  }


  if (have_male && !have_female) {
    stop("`male_candidates` supplied but `female_candidates` is missing.", call. = FALSE)
  }
  if (!have_male && have_female) {
    stop("`female_candidates` supplied but `male_candidates` is missing.", call. = FALSE)
  }

  .validate_ids(male_candidates, "male_candidates", min_len = 1)
  .validate_ids(female_candidates, "female_candidates", min_len = 1)

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
