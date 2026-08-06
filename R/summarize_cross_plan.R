#' Summarize the expected performance of a crossing plan
#'
#' Calculates the mean predicted performance of the crosses in
#' \code{cross.plan} for each available trait and value metric, including GEBV,
#' TGV, SPV, TSPV, and OHV. Crosses are matched to \code{cross.df} using their
#' parental identifiers.
#'
#' Optionally, the mean GEBV response can be plotted. When
#' \code{reference.pop} is supplied, responses are expressed relative to the
#' reference-population means. When trait \code{weights} are also supplied, the
#' corresponding desired-response direction is shown for comparison.
#'
#' @param cross.plan Data frame, matrix, or object coercible to a data frame
#'   containing the crossing plan. It must have two columns for two-way crosses
#'   or four columns for four-way crosses.
#' @param cross.df Data frame, matrix, or object coercible to a data frame
#'   containing predictions for the candidate crosses. Trait-specific value
#'   columns must follow the naming convention \code{<metric>.<trait>}, for
#'   example \code{GEBV.Yield} or \code{SPV.Height}. Crosses are matched using
#'   columns \code{parent1} and \code{parent2}, \code{male} and \code{female},
#'   or the first two columns. Four-way crosses require columns
#'   \code{parent1} to \code{parent4}.
#' @param plot Logical. If \code{TRUE}, plot the predicted mean GEBV response.
#' @param reference.pop Optional numeric matrix or data frame with individuals
#'   in rows and traits in columns. When supplied, predicted responses are
#'   expressed relative to the corresponding trait means. It must contain at
#'   least two individuals and all traits represented in \code{cross.df}.
#' @param weights Optional numeric vector of trait weights, in the same order as
#'   the traits in the returned summary. When supplied together with
#'   \code{reference.pop}, it is used to calculate and display the desired
#'   response direction.
#'
#' @return A data frame with one row per trait. The output contains:
#' \describe{
#'   \item{\code{trait}}{Trait identifier parsed from the prediction-column
#'   names.}
#'   \item{\code{ncrosses}}{Number of crosses in \code{cross.plan}.}
#'   \item{\code{mean.GEBV}, \code{mean.TGV}, \code{mean.SPV},
#'   \code{mean.TSPV}, \code{mean.OHV}}{Mean value of each available metric
#'   across the selected crosses. Metrics not available for a particular trait
#'   are returned as \code{NA}.}
#' }
#'
#' @export

summarize_cross_plan <- function(cross.plan , cross.df,plot=T, reference.pop=NA,weights=NA) {
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
      any(
        grepl(
          paste0("^", m, "\\..+$"),
          names(cross.df),
          ignore.case = TRUE
        )
      )
    }, logical(1))
  ]

  if (!length(metrics_present)) {
    .stopf(
      "No value metrics found among %s.",
      paste(all_metrics, collapse = ", ")
    )
  }

  # --- trait discovery from columns such as GEBV.name1 or SPV.name2
  cn <- names(cross.df)

  rx <- paste0(
    "^(",
    paste(metrics_present, collapse = "|"),
    ")\\.(.+)$"
  )

  mm <- regmatches(cn, regexec(rx, cn, ignore.case = TRUE))
  mm <- mm[lengths(mm) > 0]

  if (!length(mm)) {
    .stopf(
      "No trait-specific value columns found; expected names such as %s.name1.",
      paste(metrics_present, collapse = "/")
    )
  }

  # The third regex element is the part following the metric and dot
  trait_ids <- sort(unique(
    vapply(mm, `[[`, character(1), 3)
  ))

  if (!length(trait_ids)) {
    .stopf("Could not parse any trait names from value columns.")
  }

  n_crosses <- nrow(crosses)

  # helper: case-insensitive column lookup
  .find_col <- function(name) {
    hit <- names(cross.df)[tolower(names(cross.df)) == tolower(name)]

    if (length(hit)) hit[1] else NA_character_
  }

  res <- lapply(trait_ids, function(ti) {
    row <- list(
      trait = ti,
      ncrosses = n_crosses
    )

    for (m in metrics_present) {
      want <- paste0(m, ".", ti)
      hit  <- .find_col(want)

      row[[paste0("mean.", m)]] <- if (is.na(hit)) {
        NA_real_
      } else {
        v <- suppressWarnings(as.numeric(cross.df[[hit]]))

        if (all(is.na(v))) {
          NA_real_
        } else {
          mean(v, na.rm = TRUE)
        }
      }
    }

    as.data.frame(
      row,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, res)
  rownames(out) <- NULL

  if (plot) {

    # ---- Basic checks ------------------------------------------------------
    if (!"mean.GEBV" %in% names(out)) {
      .stopf(
        "Plotting requires a `mean.GEBV` column, but no GEBV columns were found."
      )
    }

    if (anyNA(out$trait) || anyDuplicated(out$trait)) {
      .stopf("Trait names must be non-missing and unique for plotting.")
    }

    if (!is.numeric(out$mean.GEBV)) {
      .stopf("`mean.GEBV` must be numeric.")
    }

    temp <- out[, c("trait", "mean.GEBV"), drop = FALSE]

    reference_missing <-
      is.null(reference.pop) ||
      (length(reference.pop) == 1L && is.atomic(reference.pop) &&
         is.na(reference.pop))

    weights_missing <-
      is.null(weights) ||
      (length(weights) == 1L && is.atomic(weights) && is.na(weights)) ||
      all(is.na(weights))

    if (reference_missing) {

      message(
        "No reference population provided; plotting absolute mean GEBVs."
      )

      predicted <- temp
      predicted$point_type <- "Predicted"

      plot_out <- ggplot2::ggplot(
        predicted,
        ggplot2::aes(x = trait, y = mean.GEBV)
      ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            xend = trait,
            y = 0,
            yend = mean.GEBV
          ),
          na.rm = TRUE
        ) +
        ggplot2::geom_point(
          ggplot2::aes(shape = point_type),
          size = 2.5,
          na.rm = TRUE
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_shape_manual(
          values = c("Predicted" = 16)
        ) +
        ggplot2::theme_grey(base_size = 10) +
        ggplot2::labs(
          x = "Traits",
          y = "Predicted population GEBV",
          title = "Predicted population GEBV",
          shape = NULL
        )

    } else {

       if (!is.data.frame(reference.pop) && !is.matrix(reference.pop)) {
        .stopf("`reference.pop` must be a numeric data frame or matrix.")
      }

      reference.pop <- as.data.frame(reference.pop)

      non_numeric <- !vapply(reference.pop, is.numeric, logical(1))

      if (any(non_numeric)) {
        .stopf(
          "All columns of `reference.pop` must be numeric. Non-numeric columns: %s.",
          paste(names(reference.pop)[non_numeric], collapse = ", ")
        )
      }

      if (nrow(reference.pop) < 2L) {
        .stopf(
          "`reference.pop` must contain at least two individuals."
        )
      }

      missing_traits <- setdiff(temp$trait, names(reference.pop))

      if (length(missing_traits)) {
        .stopf(
          "Traits  missing from `reference.pop`: %s.",
          paste(missing_traits, collapse = ", ")
        )
      }

      reference_aligned <- reference.pop[, temp$trait, drop = FALSE]

      reference_points <- colMeans(
        reference_aligned,
        na.rm = TRUE
      )

      if (any(!is.finite(reference_points))) {
        .stopf(
          "Could not calculate finite reference means for all traits."
        )
      }

      temp$mean.GEBV <- temp$mean.GEBV -
        unname(reference_points[temp$trait])

      predicted <- temp
      predicted$point_type <- "Predicted"
       if (weights_missing) {

        message(
          "No weights provided; plotting predicted responses without objective points."
        )

        plot_out <- ggplot2::ggplot(
          predicted,
          ggplot2::aes(x = trait, y = mean.GEBV)
        ) +
          ggplot2::geom_segment(
            ggplot2::aes(
              xend = trait,
              y = 0,
              yend = mean.GEBV
            ),
            na.rm = TRUE
          ) +
          ggplot2::geom_point(
            ggplot2::aes(shape = point_type),
            size = 2.5,
            na.rm = TRUE
          ) +
          ggplot2::coord_flip() +
          ggplot2::scale_shape_manual(
            values = c("Predicted" = 16)
          ) +
          ggplot2::theme_grey(base_size = 10) +
          ggplot2::labs(
            x = "Traits",
            y = "Predicted response",
            title = "Predicted response",
            shape = NULL
          )

      } else {
         weights <- as.numeric(weights)

        if (length(weights) != nrow(temp)) {
          .stopf(
            "`weights` has length %d, but %d traits are being plotted.",
            length(weights),
            nrow(temp)
          )
        }

        if (any(!is.finite(weights))) {
          .stopf("All values in `weights` must be finite.")
        }
         trait_cov <- stats::var(
          reference_aligned,
          na.rm = TRUE
        )

        if (any(!is.finite(trait_cov))) {
          .stopf(
            "The trait covariance matrix contains non-finite values."
          )
        }

        d <- drop(trait_cov %*% weights)

        if (length(d) != nrow(temp)) {
          .stopf(
            "The desired-response vector does not match the number of traits."
          )
        }

        valid <- is.finite(d) & is.finite(temp$mean.GEBV)

        if (!any(valid)) {
          .stopf(
            "No finite trait responses are available for calculating the objective response."
          )
        }


        scale_factor <- drop(
          crossprod(
            d[valid],
            temp$mean.GEBV[valid]
          ) /
            crossprod(d[valid])
        )

        objective_response <- scale_factor * d

        objective <- predicted
        objective$mean.GEBV <- objective_response
        objective$point_type <- "Objective"

        plot_data <- rbind(predicted, objective)

        plot_out <- ggplot2::ggplot(
          plot_data,
          ggplot2::aes(x = trait, y = mean.GEBV)
        ) +
           ggplot2::geom_segment(
            data = predicted,
            ggplot2::aes(
              x = trait,
              xend = trait,
              y = 0,
              yend = mean.GEBV
            ),
            inherit.aes = FALSE,
            na.rm = TRUE
          ) +
          ggplot2::geom_point(
            ggplot2::aes(shape = point_type),
            size = 2.5,
            na.rm = TRUE
          ) +
          ggplot2::coord_flip() +
          ggplot2::scale_shape_manual(
            values = c(
              "Predicted" = 16,
              "Objective" = 1
            )
          ) +
          ggplot2::theme_grey(base_size = 10) +
          ggplot2::labs(
            x = "Traits",
            y = "Predicted response",
            title = "Predicted response versus desired response",
            shape = NULL
          )
      }
    }

    print(plot_out)
  }


  out

}
