#' Evaluate a cross plan by averaging value metrics across all crosses (per trait)
#'
#' Aligns \code{cross.df} to \code{crosses} (2-way or 4-way) and computes per-trait
#' means across all crosses for available metrics (GEBV, TGV, SPV, TSPV).
#'
#' @param cross.plan  data.frame with 2 columns (2-way) or 4 columns (4-way).
#' @param cross.df data.frame (or coercible to data.frame) containing per-cross
#'   trait columns such as GEBV#, TGV#, SPV#, TSPV#. If parent columns are present
#'   (parent1..parent4 or male/female), they are used for alignment.
#'
#' @return A data.frame with one row per trait and columns:
#'   \code{trait}, \code{mean_GEBV}, \code{mean_TGV}, \code{mean_SPV},
#'   \code{mean_TSPV}, \code{mean_OHV} \code{n_crosses}.
#' @export
summarize_cross_plan <- function(cross.plan , cross.df) {
  crosses <- cross.plan
  .stopf <- function(...) stop(sprintf(...), call. = FALSE)
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

  if (!is.data.frame(crosses))  crosses  <- as.data.frame(crosses)
  if (!is.data.frame(cross.df)) cross.df <- as.data.frame(cross.df)

  if (!(ncol(crosses) %in% c(2L, 4L))) .stopf("`crosses` must have 2 or 4 columns.")

  # --- align cross.df to crosses
  if (ncol(crosses) == 2L) {
    key_crosses <- paste(crosses[[1]], crosses[[2]])
    if (all(c("parent1","parent2") %in% names(cross.df))) {
      key_df <- paste(cross.df$parent1, cross.df$parent2)
    } else if (all(c("male","female") %in% names(cross.df))) {
      key_df <- paste(cross.df$male, cross.df$female)
    } else {
      key_df <- paste(cross.df[[1]], cross.df[[2]])
    }
    idx <- match(key_crosses, key_df)

  } else { # 4-way
    need <- c("parent1","parent2","parent3","parent4")
    if (!all(need %in% names(cross.df))) {
      .stopf("For 4-way crosses, `cross.df` must contain columns: %s.", paste(need, collapse = ", "))
    }
    key_crosses <- paste(crosses[[1]], crosses[[2]], crosses[[3]], crosses[[4]])
    key_df <- paste(cross.df$parent1, cross.df$parent2, cross.df$parent3, cross.df$parent4)
    idx <- match(key_crosses, key_df)
  }

  if (anyNA(idx)) .stopf("Some crosses in `crosses` were not found in `cross.df`.")
  cross.df <- cross.df[idx, , drop = FALSE]

  # --- metrics: include only those that exist anywhere in cross.df
  all_metrics <- c("GEBV", "TGV", "SPV", "TSPV", "OHV")
  metrics_present <- all_metrics[
    vapply(all_metrics, function(m) {
      any(grepl(paste0("^", m, "[0-9]+$"), names(cross.df), ignore.case = TRUE))
    }, logical(1))
  ]
  if (!length(metrics_present)) {
    .stopf("No value metrics found among %s.", paste(all_metrics, collapse = ", "))
  }

  # --- trait discovery from present metric columns
  cn <- names(cross.df)
  rx <- paste0("^(", paste(metrics_present, collapse = "|"), ")([0-9]+)$")
  mm <- regmatches(cn, regexec(rx, cn, ignore.case = TRUE))
  mm <- mm[lengths(mm) > 0]
  if (!length(mm)) .stopf("No trait-indexed value columns found (expected %s#).",
                          paste(metrics_present, collapse = "/"))

  trait_ids <- sort(unique(as.integer(vapply(mm, `[[`, character(1), 3))))
  trait_ids <- trait_ids[is.finite(trait_ids)]
  if (!length(trait_ids)) .stopf("Could not parse any trait indices from value columns.")

  n_crosses <- nrow(crosses)

  # helper: case-insensitive column lookup
  .find_col <- function(name) {
    names(cross.df)[tolower(names(cross.df)) == tolower(name)][1] %||% NA_character_
  }

  res <- lapply(trait_ids, function(ti) {
    row <- list(trait = ti, ncrosses = n_crosses)

    for (m in metrics_present) {
      want <- paste0(m, ti)
      hit  <- .find_col(want)

      row[[paste0("mean.", m)]] <- if (is.na(hit)) {
        NA_real_
      } else {
        v <- suppressWarnings(as.numeric(cross.df[[hit]]))
        mean(v, na.rm = TRUE)
      }
    }

    as.data.frame(row, check.names = FALSE)
  })

  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}
