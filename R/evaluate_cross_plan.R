#' Evaluate a cross plan by averaging value metrics across all crosses
#'
#' @param cross_plan A data.frame with two columns as produced by `create_cross_plan()`.
#' @param cross_param Either a data.frame (when covariance = FALSE) or a list with
#'   element `$cross_values` (data.frame) when covariance = TRUE, produced by calculate_*.
#'
#' @return A single-row data.frame containing the mean across all crosses for each
#'   available value metric per trait (among EGBV, ETGV, SPV, SPTV, TSPV), with
#'   columns named `mean_<METRIC><trait_index>`, plus `n_crosses`.
#' @export
evaluate_cross_plan <- function(cross_plan, cross_param) {
  
  .stopf <- function(...) stop(sprintf(...), call. = FALSE)
  
  .validate_cross_plan <- function(x) {
    if (!is.data.frame(x)) .stopf("`cross_plan` must be a data.frame.")
    if (ncol(x) != 2L)     .stopf("`cross_plan` must have exactly 2 columns.")
    cn <- tolower(names(x))
    if (!all(cn %in% c("parent1","parent2","male","female"))) {
      .stopf("`cross_plan` columns must be either (parent1,parent2) or (male,female).")
    }
    if (anyNA(x)) .stopf("`cross_plan` contains NA.")
    invisible(NULL)
  }
  
  .extract_cross_values_df <- function(obj) {
    if (is.list(obj) && !is.null(obj$cross_values)) {
      cv <- obj$cross_values
      if (!is.data.frame(cv)) .stopf("`cross_param$cross_values` must be a data.frame.")
      return(cv)
    }
    if (is.data.frame(obj)) return(obj)
    .stopf("`cross_param` must be a data.frame or a list with `$cross_values`.")
  }
  
  .select_value_cols <- function(cn) {
    # keep value metrics (case-insensitive): EGBV, ETGV, SPV, SPTV, TSPV with trailing trait index
    keep_pat  <- "^(EGBV|ETGV|SPV|SPTV|TSPV)\\d+$"
    keep <- cn[grepl(keep_pat, cn, ignore.case = TRUE)]
    # defensively drop anything that might sneak in
    drop <- Reduce(`|`, list(
      grepl("^var",  keep, ignore.case = TRUE),
      grepl("cov",   keep, ignore.case = TRUE),
      grepl("dg",    keep, ignore.case = TRUE),
      grepl("idg",   keep, ignore.case = TRUE)
    ))
    keep[!drop]
  }
  
  # --- validate inputs
  .validate_cross_plan(cross_plan)
  values_df <- .extract_cross_values_df(cross_param)
  if (!nrow(values_df)) .stopf("`cross_param` has zero rows.")
  if (nrow(values_df) != nrow(cross_plan)) {
    .stopf("Row count mismatch: nrow(cross_plan) = %d, nrow(cross_param) = %d. They must align in order.",
           nrow(cross_plan), nrow(values_df))
  }
  
  # numeric-only
  num_cols <- vapply(values_df, is.numeric, logical(1))
  vals <- values_df[, num_cols, drop = FALSE]
  if (!ncol(vals)) .stopf("`cross_param` has no numeric columns.")
  
  # select value columns to average
  keep_cols <- .select_value_cols(names(vals))
  if (!length(keep_cols)) {
    .stopf("No value columns found among {EGBV, ETGV, SPV, SPTV, TSPV}.")
  }
  
  # compute means across all crosses
  means <- colMeans(vals[, keep_cols, drop = FALSE], na.rm = TRUE)
  
  # build single-row data.frame with mean_ prefix
  out <- as.data.frame(as.list(means), check.names = FALSE)
  names(out) <- paste0("mean_", names(out))
  out$n_crosses <- nrow(cross_plan)
  out
}
