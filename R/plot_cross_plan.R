#' Plot cross plan ridge distributions per trait
#'
#' Creates ridgeline density plots of simulated cross outcomes, separated by trait.
#' For each cross and trait, values are sampled from a Normal distribution whose mean
#' is taken from the corresponding predicted value column and whose variance is taken
#' from the available variance components.
#'
#' Variance / mean handling (per trait):
#' \itemize{
#'   \item If \code{var.A#} is available, an \emph{additive} ridge is drawn using
#'         \code{mean = EGEBV#} and \code{sd = sqrt(var.A#)}.
#'   \item If \code{var.A#} and \code{var.D#} and \code{ETGV#} are available, an additional
#'         \emph{additive+dominance} ridge is drawn using.
#'   \item If only \code{var#} is available (and no \code{var.A#}), a single ridge is drawn using
#'         \code{mean = EGEBV#} and \code{sd = sqrt(var#)}.
#' }
#'
#' Point overlays are added when the corresponding columns are present:
#' \itemize{
#'   \item \code{EGEBV#} (always used when available) and \code{SPV#} are shown as points and
#'         share the additive colour family.
#'   \item \code{ETGV#} and \code{TSPV#} are shown as points and share the additive+dominance
#'         colour family.
#'   \item \code{OHV#} is shown as a black cross (\code{shape = 4}).
#' }
#'
#' A dotted vertical line indicates the per-trait average gain (mean of \code{EGEBV#} across crosses).
#' The plot is faceted by trait (up to 3 columns).
#'
#' @param crosses data.frame with 2 columns (2-way crosses) or 4 columns (4-way crosses),
#'   defining the parents in each cross. Rows are used to label and order crosses in the plot.
#' @param cross.df data.frame (or object coercible to a data.frame) as created from 
#'   the get_variance or get_optimal_haploid_value function containing per-cross
#'   trait columns such as \code{EGEBV#}, \code{ETGV#}, \code{SPV#}, \code{TSPV#}, \code{OHV#},
#'   and variance components \code{var#}, \code{var.A#}, \code{var.D#}. If \code{cross.df}
#'   contains parent columns (e.g., \code{parent1..parent4} or \code{male,female}), they are
#'   used to align rows to \code{crosses}.
#' @param traits integer vector of trait indices to plot (e.g., \code{c(1,3)}), or \code{NULL}
#'   to plot all traits found in \code{cross.df}.
#' @param nsamples integer. Number of Monte Carlo samples per cross and trait.
#' @return A \code{ggplot} object.
#' @export
#' 
plot_cross_plan <- function(crosses, cross.df, traits = NULL,
                            nsamples = 1000L) {
  .stopf <- function(...) stop(sprintf(...), call. = FALSE)
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x
  
  if (!is.data.frame(crosses)) crosses <- as.data.frame(crosses)
  if (!is.data.frame(cross.df)) cross.df <- as.data.frame(cross.df)
  
  if (!requireNamespace("ggridges", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    .stopf("Packages `ggridges` and `ggplot2` are required.")
  }
  
  .make_cross_labels <- function(df) apply(df, 1, function(r) paste(r, collapse = "\u00D7"))
  
  # --- reorder cross.df to match crosses
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
    
  } else if (ncol(crosses) == 4L) {
    need <- c("parent1","parent2","parent3","parent4")
    if (!all(need %in% names(cross.df))) {
      .stopf("For 4-way crosses, `cross.df` must contain columns: %s.",
             paste(need, collapse = ", "))
    }
    key_crosses <- paste(crosses[[1]], crosses[[2]], crosses[[3]], crosses[[4]])
    key_df      <- paste(cross.df$parent1, cross.df$parent2, cross.df$parent3, cross.df$parent4)
    idx <- match(key_crosses, key_df)
    
  } else {
    .stopf("`crosses` must have 2 or 4 columns.")
  }
  
  if (anyNA(idx)) .stopf("Some crosses in `crosses` were not found in `cross.df`.")
  cross.df <- cross.df[idx, , drop = FALSE]
  
  # --- trait discovery
  cn <- names(cross.df)
  cn2 <- cn[!grepl("^parent[1-4]$|^male$|^female$", cn, ignore.case = TRUE)]
  rx <- "^(EGEBV|ETGV|SPV|TSPV|OHV|VAR|VAR\\.[A-Za-z]+|var|var\\.[A-Za-z]+|var[A-Za-z]+)([0-9]+)$"
  mm <- regmatches(cn2, regexec(rx, cn2, ignore.case = TRUE))
  mm <- mm[lengths(mm) > 0]
  if (!length(mm)) .stopf("No trait-indexed columns found (e.g., EGEBV1, var.A1, SPV1, etc.).")
  
  trait_chr <- vapply(mm, `[[`, character(1), 3)
  trait_ids_all <- sort(unique(as.integer(trait_chr)))
  trait_ids_all <- trait_ids_all[is.finite(trait_ids_all)]
  
  if (is.null(traits)) {
    trait_ids <- trait_ids_all
  } else {
    traits <- as.integer(traits)
    if (anyNA(traits) || any(traits < 1)) .stopf("`traits` must be positive integers.")
    missing <- setdiff(traits, trait_ids_all)
    if (length(missing)) {
      .stopf("Requested traits not available: %s. Available: %s",
             paste(missing, collapse = ", "),
             paste(trait_ids_all, collapse = ", "))
    }
    trait_ids <- traits
  }
  
  # --- flexible column pickers
  .pick <- function(prefix, ti) {
    nm <- paste0(prefix, ti)
    nm[nm %in% names(cross.df)][1] %||% NA_character_
  }
  .pick_varA <- function(ti) {
    cands <- c(paste0("var.A", ti), paste0("varA", ti), paste0("VARA", ti), paste0("VAR.A", ti))
    cands[cands %in% names(cross.df)][1] %||% NA_character_
  }
  .pick_varD <- function(ti) {
    cands <- c(paste0("var.D", ti), paste0("varD", ti), paste0("VARD", ti), paste0("VAR.D", ti))
    cands[cands %in% names(cross.df)][1] %||% NA_character_
  }
  .pick_var <- function(ti) {
    cands <- c(paste0("var", ti), paste0("VAR", ti))
    cands[cands %in% names(cross.df)][1] %||% NA_character_
  }
  .as_num <- function(x) suppressWarnings(as.numeric(x))
  .safe_sd <- function(v) {
    v[!is.finite(v) | v < 0] <- 0
    sqrt(v)
  }
  .rnorm_safe <- function(n, mean, sd) {
    if (!is.finite(mean)) mean <- 0
    if (!is.finite(sd) || sd < 0) sd <- 0
    stats::rnorm(n, mean = mean, sd = sd)
  }
  
  cross_labels <- .make_cross_labels(crosses)
  n_cross <- nrow(crosses)
  
  ridge_list <- list()
  pts_list   <- list()
  vlines     <- list()
  
  for (ti in trait_ids) {
    trait_lab <- paste0("Trait ", ti)
    
    col_E <- .pick("EGEBV", ti)
    if (is.na(col_E)) next
    mu_E <- .as_num(cross.df[[col_E]])
    
    col_ET <- .pick("ETGV", ti)
    mu_ET <- if (!is.na(col_ET)) .as_num(cross.df[[col_ET]]) else rep(NA_real_, n_cross)
    
    # variance columns
    colA <- .pick_varA(ti)
    colD <- .pick_varD(ti)
    colV <- .pick_var(ti)
    
    # --- Ridge: Additive (A) : mean EGEBV, sd sqrt(var.A)
    if (!is.na(colA)) {
      sdA <- .safe_sd(.as_num(cross.df[[colA]]))
      valsA <- lapply(seq_len(n_cross), \(i) .rnorm_safe(nsamples, mu_E[i], sdA[i]))
      ridge_list[[paste(ti, "A", sep = "_")]] <- data.frame(
        trait = trait_lab,
        group = "A",
        ridge_type = "Additive",
        cross = factor(rep(cross_labels, each = nsamples), levels = rev(unique(cross_labels))),
        value = unlist(valsA, use.names = FALSE)
      )
    }
    
    # --- Ridge: A+D : mean ETGV, sd sqrt(var.A + var.D) (only if ETGV exists)
    if (!is.na(colA) && !is.na(colD) && !is.na(col_ET)) {
      vAD <- .as_num(cross.df[[colA]]) + .as_num(cross.df[[colD]])
      sdAD <- .safe_sd(vAD)
      valsAD <- lapply(seq_len(n_cross), \(i) .rnorm_safe(nsamples, mu_ET[i], sdAD[i]))
      ridge_list[[paste(ti, "AD", sep = "_")]] <- data.frame(
        trait = trait_lab,
        group = "AD",
        ridge_type = "Additive+Dominance",
        cross = factor(rep(cross_labels, each = nsamples), levels = rev(unique(cross_labels))),
        value = unlist(valsAD, use.names = FALSE)
      )
    }
    
    # --- Ridge: Total (var) : mean EGEBV, sd sqrt(var)
    # Use only if no var.A exists (or you still want it even if var.A exists — your call)
    if (is.na(colA) && !is.na(colV)) {
      sdT <- .safe_sd(.as_num(cross.df[[colV]]))
      valsT <- lapply(seq_len(n_cross), \(i) .rnorm_safe(nsamples, mu_E[i], sdT[i]))
      ridge_list[[paste(ti, "T", sep = "_")]] <- data.frame(
        trait = trait_lab,
        group = "A",              # treat as "A" color family
        ridge_type = "Total",
        cross = factor(rep(cross_labels, each = nsamples), levels = rev(unique(cross_labels))),
        value = unlist(valsT, use.names = FALSE)
      )
    }
    
    # --- Points (grouped colors)
    add_point <- function(series, colname, group) {
      pts_list[[paste(ti, series, sep = "_")]] <<- data.frame(
        trait = trait_lab,
        cross = cross_labels,
        value = .as_num(cross.df[[colname]]),
        series = series,
        group = group
      )
    }
    
    # A family: EGEBV + SPV
    add_point("EGEBV", col_E, "A")
    col_SPV <- .pick("SPV", ti)
    if (!is.na(col_SPV)) add_point("SPV", col_SPV, "A")
    
    # AD family: ETGV + TSPV
    if (!is.na(col_ET)) add_point("ETGV", col_ET, "AD")
    col_TSPV <- .pick("TSPV", ti)
    if (!is.na(col_TSPV)) add_point("TSPV", col_TSPV, "AD")
    
    # OHV: black cross, no group coloring
    col_OHV <- .pick("OHV", ti)
    if (!is.na(col_OHV)) {
      pts_list[[paste(ti, "OHV", sep = "_")]] <- data.frame(
        trait = trait_lab,
        cross = cross_labels,
        value = .as_num(cross.df[[col_OHV]]),
        series = "OHV",
        group = NA_character_
      )
    }
    
    # vertical line: average gain (based on EGEBV means)
    vlines[[as.character(ti)]] <- data.frame(
      trait = trait_lab,
      xint = mean(mu_E, na.rm = TRUE),
      llabel = "Average gain"
    )
  }
  
  sim_df <- if (length(ridge_list)) do.call(rbind, ridge_list) else NULL
  pts_df <- if (length(pts_list))  do.call(rbind, pts_list)  else NULL
  avg_df <- if (length(vlines))    do.call(rbind, vlines)    else NULL
  if (is.null(sim_df) || !nrow(sim_df)) .stopf("No ridges could be generated for the requested traits.")
  
  # Shapes: OHV is a black cross
  shape_vals <- c(EGEBV = 16, SPV = 18, ETGV = 17, TSPV = 15, OHV = 4)
  
  p <- ggplot2::ggplot(sim_df, ggplot2::aes(x = value, y = cross)) +
    { if (!is.null(avg_df))
      ggplot2::geom_vline(
        data = avg_df,
        ggplot2::aes(xintercept = xint, linetype = llabel),
        colour = "grey40", linewidth = 0.8
      ) else ggplot2::geom_blank()
    } +
    ggridges::geom_density_ridges(
      ggplot2::aes(colour = group),
      alpha = 0, scale = 0.8, rel_min_height = 1e-7, linewidth = 1.2,
      show.legend = TRUE
    ) +
    { if (!is.null(pts_df))
      ggplot2::geom_point(
        data = pts_df[pts_df$series != "OHV", , drop = FALSE],
        ggplot2::aes(x = value, y = cross, shape = series, colour = group),
        size = 2.5, alpha = 0.9, show.legend = TRUE
      ) else ggplot2::geom_blank()
    } +
    { if (!is.null(pts_df) && any(pts_df$series == "OHV"))
      ggplot2::geom_point(
        data = pts_df[pts_df$series == "OHV", , drop = FALSE],
        ggplot2::aes(x = value, y = cross, shape = series),
        colour = "black", size = 2.8, alpha = 0.95, show.legend = TRUE
      ) else ggplot2::geom_blank()
    } +
    ggplot2::facet_wrap(~ trait, ncol = min(3, length(unique(sim_df$trait))), scales = "free_x") +
    ggplot2::xlab("Cross value") +
    ggplot2::ylab("Cross combination") +
    ggplot2::theme_grey(base_size = 10) +
    ggplot2::scale_shape_manual(name = NULL, values = shape_vals, drop = TRUE) +
    ggplot2::scale_linetype_manual(name = NULL, values = c("Average gain" = "dotted")) +
    ggplot2::guides(
      colour   = ggplot2::guide_legend(order = 1, title = NULL),
      shape    = ggplot2::guide_legend(order = 2, title = NULL),
      linetype = ggplot2::guide_legend(order = 3, title = NULL)
    ) +
    ggplot2::theme(
      legend.position   = "bottom",
      legend.key        = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background  = ggplot2::element_rect(colour = "black", fill = "grey93", linewidth = 1.3),
      axis.title        = ggplot2::element_text(size = 11),
      axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
      axis.ticks        = ggplot2::element_line(),
      axis.text         = ggplot2::element_text(size = 10)
    )
  
  p
}
