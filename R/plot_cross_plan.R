#' Plot cross plan ridge distributions per trait
#'
#' Simulates per-cross distributions using EGBV as mean and varA as variance
#' (ignoring dominance variance), draws ridgeline densities, overlays EGBV (blue circles)
#' and SPV (red diamonds, if present), and adds a dotted vertical line at the per-trait
#' average EGBV. One panel per trait (3 columns), each with its own x-axis scale.
#'
#' @param cross_plan data.frame: two columns from create_cross_plan()
#' @param cross_param data.frame or list with $cross_values from calculate_variances*()
#' @param traits integer vector of trait indices to plot, or NULL for all
#' @return ggplot object
#' @export
plot_cross_plan <- function(cross_plan, cross_param, traits = NULL) {
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

  # Map columns per trait for EGBV / varA / SPV
  .colmap_by_trait <- function(cn) {
    pat <- "^(EGBV|var[_]?a|SPV)(\\d+)$"
    keep <- cn[grepl(pat, cn, ignore.case = TRUE)]
    if (!length(keep)) .stopf("Missing columns like EGBV#, varA#, SPV#.")
    m <- regexec(pat, keep, ignore.case = TRUE)
    parts <- regmatches(keep, m)
    df <- data.frame(
      full = keep,
      kind = toupper(vapply(parts, function(x) x[[2]], character(1))),
      idx  = vapply(parts, function(x) x[[3]], character(1)),
      stringsAsFactors = FALSE
    )
    df$kind[df$kind %in% c("VARA","VAR_A")] <- "VARA"
    split(df, df$idx)
  }

  .make_cross_labels <- function(cp) {
    a <- cp[[1]]; b <- cp[[2]]
    paste0("G", a, "\u00D7", "G", b)
  }

  # --- validate / extract
  .validate_cross_plan(cross_plan)
  cv <- .extract_cross_values_df(cross_param)
  if (nrow(cv) != nrow(cross_plan)) {
    .stopf("Row mismatch: nrow(cross_plan) = %d, nrow(cross_param) = %d.", nrow(cross_plan), nrow(cv))
  }

  num_cols <- vapply(cv, is.numeric, logical(1))
  if (!any(num_cols)) .stopf("`cross_param` has no numeric columns.")
  cv_num <- cv[, num_cols, drop = FALSE]

  colmap_all <- .colmap_by_trait(names(cv_num))
  all_trait_ids <- names(colmap_all)

  # trait filtering
  if (!is.null(traits)) {
    if (!is.numeric(traits) || any(!is.finite(traits)) || any(traits < 1) || any(traits != as.integer(traits))) {
      .stopf("`traits` must be a vector of positive integers (e.g., c(1,3)).")
    }
    want <- as.character(as.integer(traits))
    missing <- setdiff(want, all_trait_ids)
    if (length(missing)) {
      .stopf("Requested traits not available: %s. Available traits: %s",
             paste(missing, collapse = ", "),
             paste(sort(as.integer(all_trait_ids)), collapse = ", "))
    }
    colmap <- colmap_all[want]
  } else {
    colmap <- colmap_all
  }
  trait_ids <- names(colmap)

  cross_labels <- .make_cross_labels(cross_plan)
  n_cross <- length(cross_labels)
  n_samples <- 1000L  # internal default

  sim_list <- list()
  mean_pts <- list()
  spv_pts  <- list()
  egvb_avg_lines <- list()

  for (ti in trait_ids) {
    meta <- colmap[[ti]]
    col_E <- meta$full[meta$kind == "EGBV"]
    col_V <- meta$full[meta$kind == "VARA"]
    col_S <- meta$full[meta$kind == "SPV"]

    if (length(col_E) != 1L || length(col_V) != 1L) {
      .stopf("Trait %s needs exactly one EGBV%s and one varA%s column.", ti, ti, ti)
    }

    mu  <- cv_num[[col_E]]
    va  <- cv_num[[col_V]]; va[!is.finite(va) | va < 0] <- 0
    sdv <- sqrt(va)

    # simulate per cross (for ridges)
    vals <- lapply(seq_len(n_cross), function(i) stats::rnorm(n_samples, mu[i], sdv[i]))
    sim_list[[ti]] <- data.frame(
      trait = paste0("Trait ", ti),
      cross = factor(rep(cross_labels, each = n_samples), levels = rev(unique(cross_labels))),
      value = unlist(vals, use.names = FALSE)
    )

    # ensure consistent legend order
    legend_levels <- c("EGBV", "SPV")

    mean_pts[[ti]] <- data.frame(
      trait  = paste0("Trait ", ti),
      cross  = cross_labels,
      value  = mu,
      legend = factor("EGBV", levels = legend_levels)
    )

    if (length(col_S) == 1L) {
      spv_pts[[ti]] <- data.frame(
        trait  = paste0("Trait ", ti),
        cross  = cross_labels,
        value  = cv_num[[col_S]],
        legend = factor("SPV", levels = legend_levels)
      )
    }

    egvb_avg_lines[[ti]] <- data.frame(
      trait  = paste0("Trait ", ti),
      xint   = mean(mu, na.rm = TRUE),
      llabel = "Average EGBV"
    )
  }

  sim_df    <- do.call(rbind, sim_list)
  means_df  <- do.call(rbind, mean_pts)
  spv_df    <- if (length(spv_pts)) do.call(rbind, spv_pts) else NULL
  avg_egbv  <- do.call(rbind, egvb_avg_lines)

  if (!requireNamespace("ggridges", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    .stopf("Packages `ggridges` and `ggplot2` are required.")
  }

  # build plot
  p <- ggplot2::ggplot(sim_df, ggplot2::aes(x = value, y = cross)) +
    # per-trait average EGBV line (dotted; its own legend entry)
    ggplot2::geom_vline(
      data = avg_egbv,
      ggplot2::aes(xintercept = xint, linetype = llabel),
      colour = "grey40", linewidth = 0.8
    ) +
    # ridgelines (no legend)
    ggridges::geom_density_ridges(
      alpha = 0, scale = 0.8, rel_min_height = 1e-7, linewidth = 1.2, show.legend = FALSE
    ) +
    # EGBV points
    ggplot2::geom_point(
      data = means_df,
      ggplot2::aes(x = value, y = cross, shape = legend, colour = legend, fill = legend),
      size = 2.5, alpha = 0.9
    ) +
    # SPV points (if present)
    { if (!is.null(spv_df))
      ggplot2::geom_point(
        data = spv_df,
        ggplot2::aes(x = value, y = cross, shape = legend, colour = legend, fill = legend),
        size = 2.5, alpha = 0.9
      )
      else ggplot2::geom_blank()
    } +
    ggplot2::facet_wrap(~ trait, ncol = min(3, length(unique(sim_df$trait))), scales = "free_x") +
    ggplot2::xlab("Cross value") +
    ggplot2::ylab("Cross combination") +
    ggplot2::theme_grey(base_size = 10) +
    # unified legend via shape; hide separate colour/fill guides but keep aesthetics for draw
    ggplot2::scale_shape_manual(
      name = NULL, values = c("EGBV" = 21, "SPV" = 23), drop = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = NULL, values = c("EGBV" = "blue", "SPV" = "red"), guide = "none", drop = TRUE
    ) +
    ggplot2::scale_fill_manual(
      name = NULL, values = c("EGBV" = "blue", "SPV" = "red"), guide = "none", drop = TRUE
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL, values = c("Average EGBV" = "dotted")
    ) +
    ggplot2::guides(
      # single legend for points (shape) with correct fill/colour in keys
      shape   = ggplot2::guide_legend(order = 1, override.aes = list(
        size = 3, alpha = 1,
        fill   = c("blue", "red"),
        colour = c("blue", "red")
      )),
      linetype = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::theme(
      legend.position   = "bottom",
      legend.key        = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background  = ggplot2::element_rect(colour = "black", fill = "grey93", linewidth = 1.3),
      axis.title        = ggplot2::element_text(size = 11),
      axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 6, r = 0, b = 0, l = 0)),
      axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 0)),
      axis.ticks        = ggplot2::element_line(),
      axis.text         = ggplot2::element_text(size = 10)
    )

  return(p)
}
