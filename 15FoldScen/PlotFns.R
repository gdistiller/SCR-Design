#Functions to summarise siumulations and do an array of plots
#function to calc appropriate summaries for RB, precision and coverage plots
make.summary <- function(df) {
  df %>%
    
    # ✅ STEP 1: row-wise calculations
    mutate(
      lower = estimate - 1.96 * se,
      upper = estimate + 1.96 * se,
      covered = (lower <= truth) & (upper >= truth)
    ) %>%
    
    # ✅ STEP 2: summarise
    group_by(design, param, stratum) %>%
    summarise(
      n = n(),
      mean_est = mean(estimate, na.rm = TRUE),
      sd_est = sd(estimate, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      sd_se = sd(se, na.rm = TRUE),
      truth = first(truth),
      
      # ✅ coverage (now works)
      coverage = mean(covered, na.rm = TRUE),
      cov_se = sqrt(coverage * (1 - coverage) / n),
      
      # Precision
      RSE_emp = sd_est / truth,
      RSE_model = mean_se / truth,
      
      # ✅ NEW: SE of precision measures
      RSE_emp_se   = (sd_est / truth) / sqrt(2 * (n - 1)),
      RSE_model_se = (sd_se / sqrt(n)) / truth,
      
      # Bias
      RB = (mean_est - truth) / truth,
      RB_se = sd_est / (sqrt(n) * truth),
      
      .groups = "drop"
    ) %>%
    
    # ✅ STEP 3: formatting / labels (
    mutate(
      design_group = case_when(
        grepl("Grid", design) ~ "Grid",
        grepl("Cluster", design) ~ "Cluster",
        grepl("Lacework", design) ~ "Lacework",
        grepl("GA4", design) ~ "GA4",
        grepl("GA5", design) ~ "GA5",
        grepl("Two Stage", design) ~ "2 Stage",
        TRUE ~ "Other"
      ),
      
      design_group = factor(
        design_group,
        levels = c("Grid", "Cluster", "Lacework", "GA4", "GA5", "2 Stage")
      ),
      
      param = case_when(
        param == "D"  ~ "Density",
        param == "L0" ~ "lambda[0]",
        param == "Sig" ~ "sigma"
      ),
      
      design = factor(
        design,
        levels = c(
          "Grid (OS)", "Grid 800m",
          "Cluster (OS)", "Cluster (2 Sig)",
          "Lacework", "Lacework (F)",
          "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both",
          "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both",
          "Two Stage"
        ),
        labels = c(
          "Grid (OS)", "Grid 800",
          "Cl (OS)", "Cl (2 Sig)",
          "LW", "LW (f)",
          "GA4 G1", "GA4 G2", "GA4 Avg", "GA4 Both", 
          "GA5 G1", "GA5 G2", "GA5 Avg", "GA5 Both",
          "2 Stage"
        )
      )
    )
}

#allows one to facet by trap level or plot with two symbols
#this is the current one
Metric.plot <- function(sum.df, plot.title, param_select, metric = "RB", ylims = NULL,
                        facet_traps = TRUE) {
  
  df <- sum.df %>%
    filter(param == param_select)
  
  # ✅ metric logic (unchanged)
  if (metric == "RB") {
    
    if (param_select == "Density") {
      ylab <- "Relative Bias (Density)"
    } else if (param_select == "lambda[0]") {
      ylab <- expression(Relative~Bias~"(" * lambda[0] * ")")
    } else if (param_select == "sigma") {
      ylab <- expression(Relative~Bias~"(" * sigma * ")")
    }
    
    yvar <- "RB"
    sevar <- "RB_se"
    ref_line <- 0
    
  } else if (metric == "RSE") {
    
    if (param_select == "Density") {
      ylab <- "Relative SE (Density)"
    } else if (param_select == "lambda[0]") {
      ylab <- expression(Relative~SE~"(" * lambda[0] * ")")
    } else if (param_select == "sigma") {
      ylab <- expression(Relative~SE~"(" * sigma * ")")
    }
    
    yvar <- "RSE_emp"
    sevar <- "RSE_emp_se"
    ref_line <- NULL
    
  } else if (metric == "coverage") {
    
    if (param_select == "Density") {
      ylab <- "Coverage (Density)"
    } else if (param_select == "lambda[0]") {
      ylab <- expression(Coverage~"(" * lambda[0] * ")")
    } else if (param_select == "sigma") {
      ylab <- expression(Coverage~"(" * sigma * ")")
    }
    
    yvar <- "coverage"
    sevar <- "cov_se"
    ref_line <- 0.95
  }
  
  # ✅ aesthetics depend on mode
  if (facet_traps) {
    base_aes <- aes(x = design,
                    y = .data[[yvar]],
                    colour = design_group)
  } else {
    base_aes <- aes(x = design,
                    y = .data[[yvar]],
                    colour = design_group,
                    shape = traps,
                    group = interaction(design_group, traps))
  }
  
  p <- ggplot(df, base_aes) +
    
    # ✅ shading
    # Design-type shading
    annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = "#009E73", alpha = 0.05) +   # Grid
    
    annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "#D55E00", alpha = 0.05) +   # Cluster
    
    annotate("rect", xmin = 4.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
             fill = "#E69F00", alpha = 0.05) +   # Lacework
    
    annotate("rect", xmin = 6.5, xmax = 10.5, ymin = -Inf, ymax = Inf,
             fill = "#0072B2", alpha = 0.05) +   # GA4
    
    annotate("rect", xmin = 10.5, xmax = 14.5, ymin = -Inf, ymax = Inf,
             fill = "#56B4E9", alpha = 0.05) +   # GA5
    
    annotate("rect", xmin = 14.5, xmax = 15.5, ymin = -Inf, ymax = Inf,
             fill = "#CC79A7", alpha = 0.05) +   # Two Stage
    
    # ✅ data
    geom_point(size = 2,
               position = if (facet_traps)
                 position_identity()
               else
                 position_dodge(width = 0.5)) +
    
    geom_errorbar(
      aes(ymin = .data[[yvar]] - 2 * .data[[sevar]],
          ymax = .data[[yvar]] + 2 * .data[[sevar]]),
      width = 0.2,
      position = if (facet_traps)
        position_identity()
      else
        position_dodge(width = 0.5)
    )
  
  # ✅ faceting logic
  if (facet_traps) {
    
    p <- p + facet_grid(
      traps ~  stratum,
      labeller = labeller(
        stratum = c("S1" = "Group 1",
                    "S2" = "Group 2"),
        traps = c("40" = "40 traps",
                  "120" = "120 traps")
      )
    )
    
  } else {
    
    p <- p + facet_wrap(
      ~ stratum,
      labeller = labeller(
        stratum = c("S1" = "Group 1",
                    "S2" = "Group 2")
      )
    ) +
      
      scale_shape_manual(
        name = "Traps",
        values = c("40" = 16, "120" = 17)
      )
  }
  
  # ✅ colours
  p <- p + scale_colour_manual(
    name = "Design type",
    values = c(
      "Grid" = "#009E73",      # green
      "Cluster" = "#D55E00",   # orange
      "Lacework" = "#E69F00",  # yellow-orange
      "GA4" = "#0072B2",       # blue
      "GA5" = "#56B4E9",       # sky blue
      "2 Stage" = "#CC79A7"    # reddish purple
    )
  ) +
    
    guides(
      colour = guide_legend(
        nrow = 1,
        byrow = TRUE
      )
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10, angle = 80, hjust = 1),
      legend.position = "top",
      strip.text = element_text(size = 12)
    ) +
    
    labs(
      y = ylab,
      x = "Design",
      title = plot.title
    )
  
  # ✅ reference line
  if (!is.null(ref_line)) {
    p <- p + geom_hline(yintercept = ref_line,
                        linetype = "dashed")
  }
  
  # ✅ limits
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  return(p)
}

#fn to combine plots with patchwork
#removes axis labels from top two
Combine.plots <- function(p1, p2, p3, global_title = NULL, tag_levels = "A", legend_position = "top", compact = TRUE) {

  # ---------------------------
  # 1. Apply compact theme to each plot separately
  # ---------------------------
  
  if (compact) {
    
    base_theme <- theme(
      legend.position = legend_position,
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      
      axis.text.x = element_text(size = 7, angle = 70, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 9),
      
      strip.text = element_text(size = 9),
      
      panel.spacing = unit(0.3, "lines"),
      plot.margin = margin(2, 2, 2, 2)
    )
    
    p1 <- p1 + base_theme
    p2 <- p2 + base_theme
    p3 <- p3 + base_theme
  }
  
  # ---------------------------
  # 2. Remove x-axis from top two plots ONLY
  # ---------------------------
  
  remove_x <- theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
  
  p1 <- p1 + remove_x
  p2 <- p2 + remove_x
  
  # p3 keeps its x-axis ✅
  
  # ---------------------------
  # 3. Combine plots
  # ---------------------------
  
  p_comb <- (p1 / p2 / p3) +
    plot_layout(
      guides = "collect",
      heights = c(1, 1, 1)
    )
  
  # ✅ 4. SET LEGEND POSITION HERE (critical fix)
  # ---------------------------
  
  p_comb <- p_comb &
    theme(legend.position = legend_position)
  
  # ---------------------------
  # 5. Add title and tags
  # ---------------------------
  
  p_comb <- p_comb +
    plot_annotation(
      title = global_title,
      tag_levels = tag_levels
    )
  
  return(p_comb)
}

#simpler version to combine the layout plots
Combine.layoutplots <- function(p1, p2, p3, global_title = NULL, tag_levels = "A", heights = c(1, 1, 1)) {
  
  # ---------------------------
  # 1. Apply compact theme
  # ---------------------------
  base_theme <- theme(
    strip.text = element_text(size = 9),
    
    panel.spacing = unit(0.25, "lines"),
    plot.margin = margin(2, 2, 2, 2)
  )
  
  p1 <- p1 + base_theme
  p2 <- p2 + base_theme
  p3 <- p3 + base_theme
  
  # ---------------------------
  # 2. Combine plots
  # ---------------------------
  
  p_comb <- (p1 / p2 / p3) +
    plot_layout(
      guides = "collect",
      heights = heights,
      widths = rep(1, 5)
    )
  
  # ---------------------------
  # 3. Title and tags
  # ---------------------------
  
  p_comb <- p_comb +
    plot_annotation(
      title = global_title,
      tag_levels = tag_levels
    )
  
  return(p_comb)
}

#plots of panels of realised designs
#latest function that allows one to either zoom in on the common array or use the full area
plot.design <- function(design.df, mask, buffer.prop = 0.2, view = c("crop", "full"), point.size = 1, 
                        trap.colour = "#d95f02", mask.size = 0.4, symbol = 16, ndim1 = 1, ndim2 = 4,
                        title.expr = NULL, levels_all = NULL) {
  
  view <- match.arg(view)
  
  mask_df <- as.data.frame(mask)
  
  # ---------------------------
  # 1. Full mask extent
  # ---------------------------
  full_mask_min_x <- min(mask_df$x)
  full_mask_max_x <- max(mask_df$x)
  full_mask_min_y <- min(mask_df$y)
  full_mask_max_y <- max(mask_df$y)
  
  full_mask_width  <- full_mask_max_x - full_mask_min_x
  full_mask_height <- full_mask_max_y - full_mask_min_y
  
  # ---------------------------
  # 2. Bounds
  # ---------------------------
  if (view == "full") {
    
    bounds <- design.df %>%
      distinct(design) %>%
      mutate(
        xmin = full_mask_min_x,
        xmax = full_mask_max_x,
        ymin = full_mask_min_y,
        ymax = full_mask_max_y
      )
    
  } else {
    
    global_bounds <- design.df %>%
      summarise(
        xmin = min(x),
        xmax = max(x),
        ymin = min(y),
        ymax = max(y)
      ) %>%
      mutate(
        xbuf = buffer.prop * (xmax - xmin),
        ybuf = buffer.prop * (ymax - ymin),
        xmin = xmin - xbuf,
        xmax = xmax + xbuf,
        ymin = ymin - ybuf,
        ymax = ymax + ybuf
      )
    
    bounds <- design.df %>%
      distinct(design) %>%
      crossing(global_bounds)
  }
  
  bounds <- bounds %>%
    left_join(design.df %>% distinct(design, design_label), by = "design")
  
  # ---------------------------
  # 3. Build mask_rep
  # ---------------------------
  mask_rep <- mask_df %>%
    crossing(bounds)

  if (!is.null(levels_all)) {
    # remove ALL dummy-expanded rows
    mask_rep <- mask_rep %>%
      filter(design_label != "dummy")
    
    # add exactly ONE dummy row back
    dummy_row <- data.frame(
      x = full_mask_min_x,
      y = full_mask_min_y,
      design = "dummy",
      design_label = "dummy",
      xmin = full_mask_min_x,
      xmax = full_mask_max_x,
      ymin = full_mask_min_y,
      ymax = full_mask_max_y,
      stringsAsFactors = FALSE
    )
    
    mask_rep <- bind_rows(mask_rep, dummy_row)
  }
  
  mask_plot <- mask_rep %>%
    filter(
      x >= xmin, x <= xmax,
      y >= ymin, y <= ymax
    )
  
  # ---------------------------
  # ✅ 4. LABEL TABLE (
  # ---------------------------
  if (!is.null(levels_all)) {
    
    label_df <- data.frame(
      design_label = levels_all,
      strip_label  = as.character(levels_all),
      stringsAsFactors = FALSE
    )
    
    design.df <- design.df %>% 
      filter(design_label %in% label_df$design_label) %>%
      left_join(label_df, by = "design_label")
    
    mask_rep <- mask_rep %>%
      filter(design_label %in% label_df$design_label) %>%
      left_join(label_df, by = "design_label")
    
    mask_plot <- mask_plot %>%
      filter(design_label %in% label_df$design_label) %>%
      left_join(label_df, by = "design_label")
    
    strip_levels <- label_df$strip_label
    
  } else {
    
    # No dummy — strip_label is just design_label as character
    design.df$strip_label <- design.df$design_label
    mask_rep$strip_label  <- mask_rep$design_label
    mask_plot$strip_label <- mask_plot$design_label
    
    strip_levels <- levels(design.df$design_label)
    
  }
  
  design.df$strip_label <- factor(design.df$strip_label, levels = strip_levels)
  mask_rep$strip_label  <- factor(mask_rep$strip_label,  levels = strip_levels)
  mask_plot$strip_label <- factor(mask_plot$strip_label, levels = strip_levels)

  # ---------------------------
  # 5. Plot
  # ---------------------------
  p <- ggplot() +
    
    # ✅ define facet structure
    geom_blank(
      data = mask_rep,
      aes(x, y)
    ) +
    
    geom_point(
      data = mask_plot %>% filter(design_label != "dummy"),
      aes(x, y),
      colour = "grey88",
      size = mask.size
    ) +
    
    geom_point(
      data = design.df %>% filter(design_label != "dummy"),
      aes(x, y),
      colour = trap.colour,
      shape = symbol,
      size = point.size
    )
  
  if (!is.null(levels_all)) {
    p <- p + geom_point(
      data = mask_rep %>% filter(design_label == "dummy"),
      aes(x, y),
      alpha = 0
    )
  }
  
  # ---------------------------
  # ✅ FACET (FINAL STABLE VERSION)
  # ---------------------------
  
  # Named display vector for labeller
  strip_display <- setNames(
    ifelse(strip_levels == "dummy", "''", strip_levels),
    strip_levels
  )
  
  p <- p +
    facet_grid(
      . ~ strip_label,
      scales = "fixed",
      space = "fixed",
      labeller = as_labeller(strip_display, default = label_parsed)
    ) +
    coord_equal() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank()
    )
  
  if (!is.null(title.expr)) {
    p <- p + labs(title = title.expr)
  }

  return(p)  
}

#scatterplot fn to construct scatterplots of RB or RSE
#note that colour scale and guide defined before
make_metric_scatter <- function(
    data,
    metric = c("RB", "RSE_emp"),
    trap_level,
    parameter = "Density",
    exclude_designs = NULL,
    label_size = 3
) {
  
  metric <- match.arg(metric)
  
  plot.df <- data |>
    filter(
      param == parameter,
      traps == trap_level
    )
  
  if (!is.null(exclude_designs)) {
    plot.df <- plot.df |>
      filter(!design %in% exclude_designs)
  }
  
  plot.df <- plot.df |>
    select(
      design,
      design_group,
      stratum,
      value = all_of(metric)
    ) |>
    pivot_wider(
      names_from = stratum,
      values_from = value
    ) |>
    left_join(design_key, by = "design") |>
    arrange(design_group, design)
  
  p <- ggplot(
    plot.df,
    aes(
      x = S1,
      y = S2,
      colour = design_group
    )
  ) +
    
    geom_point(size = 2) +
    
    geom_text_repel(
      aes(label = short_label),
      size = label_size,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    design_colour_scale +
    design_guide +
    
    theme_bw() +
    
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  if(metric == "RB") {
    
    p <- p +
      
      geom_hline(
        yintercept = 0,
        linetype = 2,
        colour = "grey50",
        show.legend = FALSE
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = 2,
        colour = "grey50",
        show.legend = FALSE
      ) +
      
      coord_cartesian(
        xlim = c(-0.1, 0.8),
        ylim = c(-0.1, 0.8)
      ) +
      
      labs(
        title = paste(trap_level, "traps"),
        x = "Relative Bias (S1)",
        y = "Relative Bias (S2)"
      )
    
  } else {
    
    p <- p +
      
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2,
        colour = "grey50",
        show.legend = FALSE
      ) +
      
      coord_cartesian(
        xlim = c(0, 2),
        ylim = c(0, 2)
      ) +
      
      labs(
        title = paste(trap_level, "traps"),
        x = "Relative Standard Error (S1)",
        y = "Relative Standard Error (S2)"
      )
  }
  
  p
}

#simplified fn to just do coverage plots for one trap lvl
Coverage.plot <- function(
    sum.df,
    trap_level,
    param_select = "Density",
    plot.title = NULL,
    ylims = NULL
) {
  
  df <- sum.df %>%
    filter(
      param == param_select,
      traps == trap_level
    ) %>%
    left_join(design_key, by = "design") %>%
    mutate(
      short_label = factor(
        short_label,
        levels = design_key$short_label
      )
    )

  if (param_select == "Density") {
    
    ylab <- "Coverage (Density)"
    
  } else if (param_select == "lambda[0]") {
    
    ylab <- expression(Coverage~"(" * lambda[0] * ")")
    
  } else if (param_select == "sigma") {
    
    ylab <- expression(Coverage~"(" * sigma * ")")
    
  } else {
    
    ylab <- "Coverage"
  }
  
  p <- ggplot(
    df,
    aes(
      x = short_label,
      y = coverage,
      colour = design_group
    )
  ) +
    
    geom_point(size = 2) +
    
    geom_hline(
      yintercept = 0.95,
      linetype = "dashed",
      show.legend = FALSE
    ) +
    
    facet_wrap(
      ~ stratum,
      nrow = 1,
      labeller = labeller(
        stratum = c(
          "S1" = "Group 1",
          "S2" = "Group 2"
        )
      )
    ) +
    
    scale_colour_manual(
      values = design_cols,
      guide = "none"
    ) +
    
    labs(
      title = plot.title,
      x = "Design",
      y = ylab
    ) +
    
    theme_bw() +
    
    theme(
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1,
        size = 7
      ),
      strip.text = element_text(size = 11),
      legend.position = "none"
    )
  
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  p
}

#new combine fn
Combine.performance.plots <- function(
    p1, p2,
    p3, p4,
    p5, p6,
    global_title = NULL,
    tag_levels = "A",
    heights = c(1,1,1)
) {
  
  p_comb <-
    (
      (p1 | p2) /
        (p3 | p4) /
        (p5 | p6)
    ) +
    plot_layout(
      guides = "collect",
      heights = heights
    ) +
    plot_annotation(
      title = global_title,
      tag_levels = tag_levels
    )
  
  p_comb &
    theme(
      legend.position = "bottom"
    )
}

#################################################
#older fns

#plots of panels of realised designs
#first one zooms into the trap array to reveal detail
plot.design.zoom <- function(design.df, mask, buffer.prop = 0.2, point.size = 1, trap.colour = "#d95f02",
                             mask.size = 0.4, symbol = 16, title.expr = "Trap configurations") {
  
  # ---------------------------
  # 1. Compute bounds per design
  # ---------------------------
  
  bounds <- design.df %>%
    group_by(design) %>%
    summarise(
      xmin = min(x),
      xmax = max(x),
      ymin = min(y),
      ymax = max(y),
      .groups = "drop"
    ) %>%
    mutate(
      xbuf = buffer.prop * (xmax - xmin),
      ybuf = buffer.prop * (ymax - ymin),
      xmin = xmin - xbuf,
      xmax = xmax + xbuf,
      ymin = ymin - ybuf,
      ymax = ymax + ybuf
    )
  
  # attach bounds to designs
  design.df <- design.df %>%
    left_join(bounds, by = "design")
  
  # ---------------------------
  # 2. Prepare mask
  # ---------------------------
  
  mask_df <- as.data.frame(mask)
  
  mask_df <- mask_df %>%
    crossing(bounds) %>%
    left_join(
      design.df %>% distinct(design, design_label),
      by = "design"
    )
  
  # crop mask per design
  mask_zoom <- mask_df %>%
    filter(
      x >= xmin, x <= xmax,
      y >= ymin, y <= ymax
    )
  
  # ---------------------------
  # 3. Plot
  # ---------------------------
  
  p <- ggplot() +
    
    geom_point(
      data = mask_zoom,
      aes(x, y),
      colour = "grey88",
      size = mask.size
    ) +
    
    geom_point(
      data = design.df,
      aes(x, y),
      shape = symbol,
      colour = trap.colour,
      size = point.size
    ) +
    
    facet_wrap(~ design_label,
               scales = "free",
               ncol = 2,
               labeller = label_parsed) +
    
    theme_bw() +
    
    theme(
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 12),
      axis.ticks = element_line(),
      panel.grid = element_blank()
    ) +
    
    labs(
      x = "x coordinate",
      y = "y coordinate"
    )
  # ---------------------------
  # optional title
  # ---------------------------
  
  if (!is.null(title.expr)) {
    p <- p + labs(title = title.expr)
  }
  
  return(p)
}

#2nd one keeps the same spatial extent for all plots
plot.design.fixed <- function(design.df, mask, point.size = 2, mask.size = 0.3, trap.colour = "#d95f02", symbol = 3, ncol = 2,
                              title.expr = NULL) {
  
  # convert mask once
  mask_df <- as.data.frame(mask)
  
  # ---------------------------
  # build plot
  # ---------------------------
  
  p <- ggplot() +
    
    # mask
    geom_point(
      data = mask_df,
      aes(x, y),
      colour = "grey90",
      size = mask.size
    ) +
    
    # traps
    geom_point(
      data = design.df,
      aes(x, y),
      colour = trap.colour,
      shape = symbol,
      size = point.size
    ) +
    
    # facets (fixed extent)
    facet_wrap(~ design_label,
               scales = "fixed",
               labeller = label_parsed,
               dir = "h",
               ncol = ncol) +
    
    # spatial accuracy
    coord_equal() +
    
    # theme
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 9),
      strip.text = element_text(size = 12),
      panel.grid = element_blank()
    ) +
    
    labs(
      x = "x",
      y = "y"
    )
  
  # ---------------------------
  # optional title
  # ---------------------------
  
  if (!is.null(title.expr)) {
    p <- p + labs(title = title.expr)
  }
  
  return(p)
}

plot.design.inset = function(mask, traps, inset_box_size = 0.2, inset_pad = 0.06, crop_buffer = 6000){
  
  ## Full mask extent, calculated from mask
  full_mask_min_x <- min(mask$x)
  full_mask_max_x <- max(mask$x)
  full_mask_min_y <- min(mask$y)
  full_mask_max_y <- max(mask$y)
  
  full_mask_width  <- full_mask_max_x - full_mask_min_x
  full_mask_height <- full_mask_max_y - full_mask_min_y
  
  
  ## Trap bounding box
  trap_min_x <- min(traps$x)
  trap_max_x <- max(traps$x)
  trap_min_y <- min(traps$y)
  trap_max_y <- max(traps$y)
  
  
  ## Cropped mask extent:
  ## trap bounding box + buffer, constrained to full mask extent
  cropped_mask_min_x <- max(full_mask_min_x, trap_min_x - crop_buffer)
  cropped_mask_max_x <- min(full_mask_max_x, trap_max_x + crop_buffer)
  cropped_mask_min_y <- max(full_mask_min_y, trap_min_y - crop_buffer)
  cropped_mask_max_y <- min(full_mask_max_y, trap_max_y + crop_buffer)
  
  crop_width  <- cropped_mask_max_x - cropped_mask_min_x
  crop_height <- cropped_mask_max_y - cropped_mask_min_y
  
  ## Crop the mask points actually shown in the zoomed plot
  cropped_mask_xy <- mask %>%
    dplyr::filter(
      x >= cropped_mask_min_x,
      x <= cropped_mask_max_x,
      y >= cropped_mask_min_y,
      y <= cropped_mask_max_y
    )
  
  ## Inset position in the zoomed plot
  inset_margin_x <- 0.04 * crop_width
  inset_margin_y <- 0.04 * crop_height
  
  inset_outer_width  <- inset_box_size * crop_width
  inset_outer_height <- inset_outer_width * full_mask_height / full_mask_width
  
  inset_outer <- data.frame(
    xmin = cropped_mask_max_x - inset_margin_x - inset_outer_width,
    xmax = cropped_mask_max_x - inset_margin_x,
    ymin = cropped_mask_max_y - inset_margin_y - inset_outer_height,
    ymax = cropped_mask_max_y - inset_margin_y
  )
  
  ## Inner grey box: full mask extent, drawn inside the inset border
  inset_inner <- data.frame(
    xmin = inset_outer$xmin + inset_pad * inset_outer_width,
    xmax = inset_outer$xmax - inset_pad * inset_outer_width,
    ymin = inset_outer$ymin + inset_pad * inset_outer_height,
    ymax = inset_outer$ymax - inset_pad * inset_outer_height
  )
  
  inset_inner_width  <- inset_inner$xmax - inset_inner$xmin
  inset_inner_height <- inset_inner$ymax - inset_inner$ymin
  
  ## Red box: cropped-mask extent, mapped into the full-mask inset
  inset_crop <- data.frame(
    xmin = inset_inner$xmin +
      ((cropped_mask_min_x - full_mask_min_x) / full_mask_width) * inset_inner_width,
    
    xmax = inset_inner$xmin +
      ((cropped_mask_max_x - full_mask_min_x) / full_mask_width) * inset_inner_width,
    
    ymin = inset_inner$ymin +
      ((cropped_mask_min_y - full_mask_min_y) / full_mask_height) * inset_inner_height,
    
    ymax = inset_inner$ymin +
      ((cropped_mask_max_y - full_mask_min_y) / full_mask_height) * inset_inner_height
  )
  
  ## Plot of full area
  ggplot() +
    geom_point(
      data = mask,
      shape = 21,
      aes(x = x, y = y),
      colour = "grey80",
      fill = "grey80",
      size = 1
    ) +
    geom_point(
      data = traps,
      aes(x = x, y = y),
      shape = 4,
      colour = "red",
      size = 2,
      stroke = 1.2
    ) +
    coord_equal(
      xlim = c(full_mask_min_x, full_mask_max_x),
      ylim = c(full_mask_min_y, full_mask_max_y),
      expand = FALSE
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank()
    )
  
  ## Plot of cropped area, with dynamic inset
  p0 <- ggplot() +
    geom_point(
      data = cropped_mask_xy,
      shape = 21,
      aes(x = x, y = y),
      colour = "grey80",
      fill = "grey80",
      size = 1
    ) +
    geom_point(
      data = traps,
      aes(x = x, y = y),
      shape = 4,
      colour = "red",
      size = 2,
      stroke = 1.2
    ) +
    
    ## inset outer border
    geom_rect(
      data = inset_outer,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = "black",
      linewidth = 0.6
    ) +
    
    ## inset full mask extent
    geom_rect(
      data = inset_inner,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey80",
      colour = NA
    ) +
    
    ## inset cropped-mask extent
    geom_rect(
      data = inset_crop,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      colour = "red",
      linewidth = 0.8
    ) +
    
    coord_equal(
      xlim = c(cropped_mask_min_x, cropped_mask_max_x),
      ylim = c(cropped_mask_min_y, cropped_mask_max_y),
      expand = FALSE
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank()
    )
  p0
}

#revised general plotting function to plot one measure at a time for a given parameter and trap level
Metric.plot <- function(sum.df, plot.title, param_select, metric = "RB", ylims = NULL) {
  
  # ✅ filter parameter a
  df <- sum.df %>%
    filter(param == param_select)
  
  # ✅ choose metric + SE + reference line
  if (metric == "RB") {
    
    if (param_select == "Density") {
      ylab <- "Relative Bias (Density)"
    } else if (param_select == "lambda[0]") {
      ylab <- expression("Relative Bias (" * lambda[0] * ")")
    } else if (param_select == "sigma") {
      ylab <- expression("Relative Bias (" * sigma * ")")
    }
    
    yvar <- "RB"
    sevar <- "RB_se"
    ref_line <- 0
    
  } else if (metric == "RSE") {
    
    if (param_select == "Density") {
      ylab <- "Relative SE (Density)"
    } else if (param_select == "lambda[0]") {
      ylab <- expression("Relative SE (" * lambda[0] * ")")
    } else if (param_select == "sigma") {
      ylab <- expression("Relative SE (" * sigma * ")")
    }
    
    yvar <- "RSE_emp"
    sevar <- "RSE_emp_se"
    ref_line <- NULL
    
  } else if (metric == "coverage") {
    
    if (param_select == "Density") {
      ylab <- "Coverage (Density)"
    } else if (param_select == "lambda[0]") {
      ylab <- expression("Coverage (" * lambda[0] * ")")
    } else if (param_select == "sigma") {
      ylab <- expression("Coverage (" * sigma * ")")
    }
    
    yvar <- "coverage"
    sevar <- "cov_se"
    ref_line <- 0.95
  }
  
  # ✅ base plot
  p <- ggplot(df,
              aes(x = design,
                  y = .data[[yvar]],
                  colour = design_group)) +
    
    # ✅ updated shading (5 blocks)
    annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = "#1b9e77", alpha = 0.08) +
    annotate("rect", xmin = 2.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
             fill = "#d95f02", alpha = 0.08) +
    annotate("rect", xmin = 5.5, xmax = 9.5, ymin = -Inf, ymax = Inf,
             fill = "#6a51a3", alpha = 0.08) +
    annotate("rect", xmin = 9.5, xmax = 13.5, ymin = -Inf, ymax = Inf,
             fill = "#9e9ac8", alpha = 0.08) +
    annotate("rect", xmin = 13.5, xmax = 14.5, ymin = -Inf, ymax = Inf,
             fill = "#bcbddc", alpha = 0.08) +
    
    # ✅ data
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = .data[[yvar]] - 2 * .data[[sevar]],
          ymax = .data[[yvar]] + 2 * .data[[sevar]]),
      width = 0.2
    ) +
    
    # ✅ FIXED FACETING
    facet_grid(
      ~ stratum + traps,
      labeller = labeller(
        stratum = c("S1" = "Stratum 1",
                    "S2" = "Stratum 2"),
        traps = c("40" = "40 traps",
                  "120" = "120 traps")
      )
    ) +
    
    # ✅ updated colours / legend labels
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "GA4" = "#6a51a3",
        "GA5" = "#9e9ac8",
        "2 Stage" = "#bcbddc"
      )
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10, angle = 80, hjust = 1),
      legend.position = "top",
      strip.text = element_text(size = 12)
    ) +
    
    labs(
      y = ylab,
      x = "Design",
      title = plot.title
    )
  
  # ✅ add reference line if needed
  if (!is.null(ref_line)) {
    p <- p + geom_hline(yintercept = ref_line,
                        linetype = "dashed")
  }
  
  # ✅ apply y limits if provided
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  return(p)
}

#plotting functions, here RB
RB.plot <- function(sum.df, plot.title, ylims){
  ggplot(sum.df, aes(x = design, y = RB, colour = design_group)) +
    
    # ✅ Shading (4 groups, with similar colours for 3a and 3b)
    
    geom_rect(
      data = data.frame(
        xmin = c(0.5, 2.5, 5.5, 10.5),
        xmax = c(2.5, 5.5, 10.5, 15.5),
        ymin = -Inf,
        ymax = Inf,
        fill_col = c("#1b9e77", "#d95f02", "#6a51a3", "#bcbddc")
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    scale_fill_identity() +
    
    
    # ✅ Data
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = RB - 2 * RB_se, ymax = RB + 2 * RB_se),
      width = 0.2
    ) +
    
    # ✅ Facets
    facet_grid(
      param ~ stratum,
      labeller = labeller(
        param = label_parsed,
        stratum = c(
          "S1" = "Stratum 1",
          "S2" = "Stratum 2"
        )
      )
    ) +
    
    # ✅ Reference line
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    # ✅ Point colours (group-based)
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "Optimised 1" = "#6a51a3",
        "Optimised 2" = "#7570b3"
      )
    ) +
    coord_cartesian(ylim = ylims) +
    
    # ✅ Theme
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "top",
      strip.text = element_text(size = 12)   # increase from default (~10)
    ) +
    
    labs(
      y = "Relative Bias",
      x = "Design",
      title = plot.title
    )
}

#next two try to combine the two trap levels, operates on a combined df
#RB.plot 2 nests trap lvl within stratum, but x axis labels are too squashed
RB.plot2 <- function(sum.df, plot.title, ylims){
  ggplot(sum.df, aes(x = design, y = RB, colour = design_group)) +
    
    # ✅ Shading (4 groups, with similar colours for 3a and 3b)
    
    geom_rect(
      data = data.frame(
        xmin = c(0.5, 2.5, 5.5, 10.5),
        xmax = c(2.5, 5.5, 10.5, 15.5),
        ymin = -Inf,
        ymax = Inf,
        fill_col = c("#1b9e77", "#d95f02", "#6a51a3", "#bcbddc")
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    scale_fill_identity() +
    
    
    # ✅ Data
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = RB - 2 * RB_se, ymax = RB + 2 * RB_se),
      width = 0.2
    ) +
    
    # ✅ Facets
    facet_grid(
      param ~ stratum + traps,
      labeller = labeller(
        param = label_parsed,
        stratum = c(
          "S1" = "Stratum 1",
          "S2" = "Stratum 2"
        ),
        traps = c("40" = "40 traps", "120" = "120 traps")
      )
    ) +
    
    # ✅ Reference line
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    # ✅ Point colours (group-based)
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "Optimised 1" = "#6a51a3",
        "Optimised 2" = "#7570b3"
      )
    ) +
    coord_cartesian(ylim = ylims) +
    
    # ✅ Theme
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "top",
      strip.text = element_text(size = 11),   # increase from default (~10),
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    
    labs(
      y = "Relative Bias",
      x = "Design",
      title = plot.title
    )
}

#this plots the two levels with different symbols, so similar to precision plot
RB.plot.combined <- function(sum.df, plot.title, ylims){
  
  ggplot(sum.df,
         aes(x = design,
             y = RB,
             colour = design_group,
             shape = traps)) +
    
    # ✅ shading stays the same
    geom_rect(
      data = data.frame(
        xmin = c(0.5, 2.5, 5.5, 10.5),
        xmax = c(2.5, 5.5, 10.5, 15.5),
        ymin = -Inf,
        ymax = Inf,
        fill_col = c("#1b9e77", "#d95f02", "#6a51a3", "#bcbddc")
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    scale_fill_identity() +
    
    # ✅ side-by-side positioning
    geom_point(size = 2,
               position = position_dodge(width = 0.5)) +
    
    geom_errorbar(
      aes(ymin = RB - 2 * RB_se,
          ymax = RB + 2 * RB_se),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    
    # ✅ facets unchanged
    facet_grid(
      param ~ stratum,
      labeller = labeller(
        param = label_parsed,
        stratum = c(
          "S1" = "Stratum 1",
          "S2" = "Stratum 2"
        )
      )
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed") +
    
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "Optimised 1" = "#6a51a3",
        "Optimised 2" = "#7570b3"
      )
    ) +
    
    # ✅ NEW: trap legend
    scale_shape_manual(
      name = "Traps",
      values = c("40" = 16, "120" = 17)
    ) +
    
    coord_cartesian(ylim = ylims) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "top",
      strip.text = element_text(size = 11),  # slightly tighter,
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    
    labs(
      y = "Relative Bias",
      x = "Design",
      title = plot.title
    )
}

#general fn to use for RB and coverage (and maybe others)
#yvar - column name as string (e.g. "RB")
#sevar: SE column (e.g. "RB_se")
#ylab: y-axis label
#yline: reference line (e.g. 0 or 0.95)
plot_metric <- function(sum.df, yvar, sevar, ylab, yline = NULL, plot.title, ylims = NULL) {
  
  ggplot(sum.df,
         aes(x = design,
             y = .data[[yvar]],
             colour = design_group,
             shape = traps,
             group = interaction(design_group, traps))) +
    
    # ✅ shading
    geom_rect(
      data = data.frame(
        xmin = c(0.5, 2.5, 5.5, 10.5),
        xmax = c(2.5, 5.5, 10.5, 15.5),
        ymin = -Inf,
        ymax = Inf,
        fill_col = c("#1b9e77", "#d95f02", "#6a51a3", "#bcbddc")
      ),
      aes(xmin = xmin, xmax = xmax,
          ymin = ymin, ymax = ymax,
          fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    scale_fill_identity() +
    
    # ✅ points
    geom_point(
      size = 2,
      position = position_dodge(width = 0.5)
    ) +
    
    # ✅ error bars
    geom_errorbar(
      aes(ymin = .data[[yvar]] - 2 * .data[[sevar]],
          ymax = .data[[yvar]] + 2 * .data[[sevar]]),
      width = 0.2,
      position = position_dodge(width = 0.5)
    ) +
    
    # ✅ optional reference line
    {if (!is.null(yline))
      geom_hline(yintercept = yline,
                 linetype = "dashed",
                 colour = "black",
                 linewidth = 0.8)
    } +
    
    # ✅ facets
    facet_grid(
      param ~ stratum,
      labeller = labeller(
        param = label_parsed,
        stratum = c("S1" = "Stratum 1",
                    "S2" = "Stratum 2")
      )
    ) +
    
    # ✅ colours
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "Optimised 1" = "#6a51a3",
        "Optimised 2" = "#7570b3"
      )
    ) +
    
    # ✅ shapes (traps)
    scale_shape_manual(
      name = "Traps",
      values = c("40" = 16, "120" = 17),
      labels = c("40" = "40 traps", "120" = "120 traps")
    ) +
    
    # ✅ limits
    {if (!is.null(ylims)) coord_cartesian(ylim = ylims)} +
    
    # ✅ theme
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
      strip.text = element_text(size = 12),
      legend.position = "top",
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    
    # ✅ labels
    labs(
      y = ylab,
      x = "Design",
      title = plot.title
    )
}

#plotting two measures of precision
Prec.plot <- function(sum.df, plot.title, ylims){
  
  #reshape so that the two Rel SEs can be compared
  precision_long <- sum.df %>%
    pivot_longer(
      cols = c(RSE_emp, RSE_model),
      names_to = "type",
      values_to = "RSE"
    )
  
  precision_long$type <- recode(precision_long$type,
                                RSE_emp = "Empirical",
                                RSE_model = "Model-based"
  )
  
  ggplot(precision_long,
         aes(x = design, y = RSE, colour = design_group, shape = type)) +
    
    # ✅ Same shading as before
    geom_rect(
      data = data.frame(
        xmin = c(0.5, 2.5, 5.5, 10.5),
        xmax = c(2.5, 5.5, 10.5, 15.5),
        ymin = -Inf,
        ymax = Inf,
        fill_col = c("#1b9e77", "#d95f02", "#6a51a3", "#bcbddc")
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    scale_fill_identity() +
    
    # ✅ Points (now two per design)
    geom_point(size = 2, position = position_dodge(width = 0.4)) +
    
    # ✅ Facets
    facet_grid(
      param ~ stratum,
      labeller = labeller(
        param = label_parsed,
        stratum = c("S1" = "Stratum 1", "S2" = "Stratum 2")
      )
    ) +
    
    # ✅ Colour scale
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "Optimised 1" = "#6a51a3",
        "Optimised 2" = "#bcbddc"
      )
    ) +
    
    # ✅ Shape scale for precision type
    scale_shape_manual(
      name = "Precision type",
      values = c(
        "Empirical" = 16,   # filled circle
        "Model-based" = 17  # triangle
      )
    ) +
    coord_cartesian(ylim = ylims) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
      strip.text = element_text(size = 12),
      legend.position = "top"
    ) +
    
    labs(
      y = "Relative Standard Error",
      x = "Design",
      title = plot.title
    )
}

#coverage plot
Coverage.plot <- function(sum.df, plot.title, ylims){
  ggplot(sum.df, aes(x = design, y = coverage, colour = design_group)) +
    
    # ✅ same shading
    geom_rect(
      data = data.frame(
        xmin = c(0.5, 2.5, 5.5, 10.5),
        xmax = c(2.5, 5.5, 10.5, 15.5),
        ymin = -Inf,
        ymax = Inf,
        fill_col = c("#1b9e77", "#d95f02", "#6a51a3", "#bcbddc")
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col),
      inherit.aes = FALSE,
      alpha = 0.08
    ) +
    scale_fill_identity() +
    
    # ✅ points
    geom_point(size = 2) +
    
    # ✅ optional error bars
    geom_errorbar(
      aes(ymin = coverage - 2 * cov_se, ymax = coverage + 2 * cov_se),
      width = 0.2
    ) +
    
    # ✅ target line
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black", linewidth = 0.8) +
    
    # ✅ facets
    facet_grid(
      param ~ stratum,
      labeller = labeller(
        param = label_parsed,
        stratum = c("S1" = "Stratum 1", "S2" = "Stratum 2")
      )
    ) +
    
    # ✅ colours
    scale_colour_manual(
      name = "Design type",
      values = c(
        "Grid" = "#1b9e77",
        "Cluster" = "#d95f02",
        "Optimised 1" = "#6a51a3",
        "Optimised 2" = "#bcbddc"
      )
    ) +
    
    coord_cartesian(ylim = ylims) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 10, angle = 60, hjust = 1),
      strip.text = element_text(size = 12),
      legend.position = "top"
    ) +
    
    labs(
      y = "Coverage probability",
      x = "Design",
      title = plot.title
    )
}

#next set of functions are for plotting the dbns of estimates in quarto
#plot_param_raw is the function called by plot_block
plot_param_raw <- function(df, design_name, param_name, se = FALSE, rel_se = TRUE, ylim = NULL) {
  
  param_labels <- c(
    D = "Density",
    L0 = "lambda[0]",
    Sig = "sigma"
  )
  
  param_expr <- function(p) {
    parse(text = param_labels[[p]])[[1]]
  }
  
  sub <- df %>%
    filter(
      design == design_name,
      param == param_name
    )
  
  # choose variable
  if (!se) {
    sub <- sub %>% mutate(y = estimate)
  } else {
    if (rel_se) {
      sub <- sub %>% mutate(y = se / truth)
    } else {
      sub <- sub %>% mutate(y = se)
    }
  }
  
  ggplot(sub, aes(x = 1, y = y)) +
    
    geom_boxplot(outlier.alpha = 0.3) +
    
    # ✅ separate panels per stratum with independent scales
    facet_wrap(~ stratum, scales = "free_y") +
    
    # ✅ truth line per panel (only for estimates)
    {if (!se) geom_hline(aes(yintercept = truth),
                         linetype = "dashed",
                         colour = "red")} +
    
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    
    labs(
      title = design_name,  # keep title simple
      
      y = if (!se) {
        bquote("Estimates of" ~ .(param_expr(param_name)))
      } else if (rel_se) {
        bquote("Relative SE of" ~ .(param_expr(param_name)))
      } else {
        bquote("SE of" ~ .(param_expr(param_name)))
      },
      x = NULL
    )
}

#cycles through the parameters for a given design
plot_block <- function(df, design_name, se = FALSE, rel_se = TRUE) {
  
  plots <- lapply(params, function(p) {
    plot_param_raw(
      df,
      design_name = design_name,
      param_name = p,
      se = se,
      rel_se = rel_se
    )
  })
  
  # ✅ Important: print plots so Quarto renders them
  for (p in plots) suppressWarnings(print(p))
}

#another version, very similar to the original, going to use original
plot_param_raw2 <- function(df, design_name, param_name, use_se = FALSE, rel_se = TRUE) {
  
  # ✅ parameter labels
  param_labels <- c(
    D = "Density",
    L0 = "lambda[0]",
    Sig = "sigma"
  )
  
  param_expr <- function(p) {
    parse(text = param_labels[[p]])[[1]]
  }
  
  # ✅ subset data
  sub <- df %>%
    filter(
      design == design_name,
      param == param_name
    )
  
  # ✅ choose variable OUTSIDE mutate (important!)
  if (!use_se) {
    sub <- sub %>% mutate(y = estimate)
  } else if (rel_se) {
    sub <- sub %>% mutate(y = se / truth)
  } else {
    sub <- sub %>% mutate(y = se)
  }
  
  # ✅ remove problematic values
  sub <- sub %>% filter(is.finite(y))
  
  # ✅ plot
  ggplot(sub, aes(x = stratum, y = y, fill = stratum)) +
    
    geom_boxplot(
      width = 0.6,
      alpha = 0.7,
      outlier.size = 0.6
    ) +
    facet_wrap(~ stratum, scales = "free_y") +
    
    # ✅ truth lines per stratum (only for estimates)
    {if (!use_se)
      geom_hline(
        aes(yintercept = truth, colour = stratum),
        linetype = "dashed",
        linewidth = 0.6,
        show.legend = FALSE
      )
    } +
    
    # ✅ nice colours
    scale_fill_manual(
      values = c("S1" = "#1b9e77", "S2" = "#d95f02"),
      labels = c("S1" = "Stratum 1", "S2" = "Stratum 2")
    ) +
    
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",         # redundant (x-axis already shows strata)
      panel.spacing = unit(0.3, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    
    labs(
      title = design_name,
      y = if (!use_se) {
        bquote("Estimates of" ~ .(param_expr(param_name)))
      } else if (rel_se) {
        bquote("Relative SE of" ~ .(param_expr(param_name)))
      } else {
        bquote("SE of" ~ .(param_expr(param_name)))
      },
      x = "Stratum"
    )
}

#plot_block fn that calls version2
plot_block2 <- function(df, design_name, use_se = FALSE, rel_se = TRUE) {
  
  plots <- lapply(params, function(p) {
    plot_param_raw2(
      df,
      design_name = design_name,
      param_name = p,
      se = use_se,
      rel_se = rel_se
    )
  })
  
  # ✅ Important: print plots so Quarto renders them
  for (p in plots) suppressWarnings(print(p))
}

#these next two operate on a df that has the two trap lvls combined
#the boxplots plot both trap lvls within each strata rather than using tabs
plot_param_combined <- function(df, design_name, se = FALSE, rel_se = TRUE) {
  
  param_labels <- c(
    D = "Density",
    L0 = "lambda[0]",
    Sig = "sigma"
  )
  
  param_expr <- function(p) {
    parse(text = param_labels[[p]])[[1]]
  }
  
  sub <- df %>%
    filter(
      design == design_name,
      param == param_name
    ) %>%
    mutate(
      y = ifelse(!se, estimate,
                 ifelse(rel_se, se / truth, se))
    ) %>%
    filter(is.finite(y))
  
  ggplot(sub, aes(x = traps, y = y, fill = traps)) +
    
    geom_boxplot(width = 0.6, outlier.size = 0.6, alpha = 0.7) +
    
    facet_wrap(~ stratum, scales = "free_y") +
    
    {if (!se)
      geom_hline(aes(yintercept = truth),
                 linetype = "dashed",
                 linewidth = 0.6,
                 colour = "black")} +
    
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 16)
    ) +
    
    labs(
      title = design_name,
      y = if (!se) {
        bquote("Estimates of " ~ .(param_expr(param_name)))
      } else if (rel_se) {
        bquote("Relative SE of " ~ .(param_expr(param_name)))
      } else {
        bquote("SE of " ~ .(param_expr(param_name)))
      },
      x = "Number of traps"
    )
}

plot_block_combined <- function(df, design_name, se = FALSE, rel_se = TRUE) {
  
  plots <- lapply(params, function(p) {
    plot_param_combined(
      df,
      design_name = design_name,
      param_name = p,
      se = se,
      rel_se = rel_se
    )
  })
  
  for (p in plots) print(p)
}

#new fn to tab the parameters rather than trap lvl
#still calls plot_param_combined to show 4 and 120 traps together
plot_block_tabs <- function(df, design_name,
                            use_se = FALSE, rel_se = TRUE) {
  
  cat("::: panel-tabset\n\n")
  
  plot_block_tabs <- function(df, design_name,
                              use_se = FALSE, rel_se = TRUE) {
    
    cat("::: panel-tabset\n\n")
    
    nice_names <- c(
      D   = "Density",
      L0  = "λ₀",
      Sig = "σ"
    )
    
    for (p in params) {
      
      # ✅ Use the nice label instead of raw "D", "L0", "Sig"
      cat("###", nice_names[[p]], "\n\n")
      
      print(
        plot_param_combined(
          df,
          design_name = design_name,
          param_name = p,
          use_se = use_se,
          rel_se = rel_se
        )
      )
      
      cat("\n\n")  # reset markdown (important for tab parsing)
    }
    
    cat(":::\n\n")
  }
  
  for (p in params) {
    
    # ✅ Use the nice label instead of raw "D", "L0", "Sig"
    cat("###", nice_names[[p]], "\n\n")
    
    print(
      plot_param_combined(
        df,
        design_name = design_name,
        param_name = p,
        use_se = use_se,
        rel_se = rel_se
      )
    )
    
    cat("\n\n")  # reset markdown (important for tab parsing)
  }
  
  cat(":::\n\n")
}

#plot to generate boxplots of estimates in a similar structure
Boxplot.general <- function(df, plot.title, yvar = "estimate", show_ref = TRUE, ref_df = NULL) {
  
  # ✅ Fix parameter labels
  df <- df %>%
    mutate(
      param_label = case_when(
        param == "D"   ~ "D",
        param == "L0"  ~ "lambda[0]",
        param == "Sig" ~ "sigma"
      )
    )
  
  # ✅ Base plot
  p <- ggplot(df,
              aes(x = design,
                  y = .data[[yvar]],
                  fill = traps,
                  group = interaction(design, traps))) +
    
    # ✅ Shading
    annotate("rect", xmin = 0.5,  xmax = 2.5,  ymin = -Inf, ymax = Inf,
             fill = "#1b9e77", alpha = 0.08) +
    annotate("rect", xmin = 2.5,  xmax = 5.5,  ymin = -Inf, ymax = Inf,
             fill = "#d95f02", alpha = 0.08) +
    annotate("rect", xmin = 5.5,  xmax = 10.5, ymin = -Inf, ymax = Inf,
             fill = "#6a51a3", alpha = 0.08) +
    annotate("rect", xmin = 10.5, xmax = 15.5, ymin = -Inf, ymax = Inf,
             fill = "#bcbddc", alpha = 0.08) +
    
    # ✅ Boxplots
    geom_boxplot(
      position = position_dodge(width = 0.75),
      alpha = 0.7,
      outlier.size = 0.5
    ) +
    
    # ✅ Facets
    facet_wrap(
      ~ param_label + stratum,
      scales = "free_y",
      ncol = 2,
      labeller = labeller(
        param_label = label_parsed,
        stratum = c(
          "S1" = "Stratum 1",
          "S2" = "Stratum 2"
        )
      )
    ) +
    
    # ✅ Trap legend
    scale_fill_manual(
      name = "Traps",
      values = c(
        "40"  = "#1b9e77",
        "120" = "#d95f02"
      ),
      labels = c(
        "40"  = "40 traps",
        "120" = "120 traps"
      )
    ) +
    
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      strip.text = element_text(size = 12),
      legend.position = "top"
    ) +
    
    labs(
      y = yvar,
      x = "Design",
      title = plot.title
    )
  
  # ✅ ✅ Add REF lines (stratum-specific)
  if (show_ref & !is.null(ref_df)) {
    
    ref_df_fixed <- ref_df %>%
      mutate(
        param_label = case_when(
          param == "D"   ~ "D",
          param == "L0"  ~ "lambda[0]",
          param == "Sig" ~ "sigma"
        )
      ) %>%
      distinct(param_label, stratum, ref_val)
    
    p <- p +
      geom_hline(
        data = ref_df_fixed,
        aes(yintercept = ref_val),
        linetype = "solid",
        colour = "grey",
        linewidth = 1,
        inherit.aes = FALSE
      )
  }
  
  return(p)
}4



