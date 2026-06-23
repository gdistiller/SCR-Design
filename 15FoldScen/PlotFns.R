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
    
    # ✅ STEP 3: formatting / labels (unchanged)
    mutate(
      design_group = case_when(
        grepl("Grid", design) ~ "Grid",
        grepl("Cluster", design) | grepl("Lacework", design) ~ "Cluster",
        grepl("GA4", design) ~ "GA4",
        grepl("GA5", design) ~ "GA5",
        grepl("Two Stage", design) ~ "GA5",
        TRUE ~ "Other"
      ),
      
      design_group = factor(
        design_group,
        levels = c("Grid", "Cluster", "GA4", "GA5", "2 Stage")
      ),
      
      param = case_when(
        param == "D"  ~ "Density",
        param == "L0" ~ "lambda[0]",
        param == "Sig" ~ "sigma"
      ),
      
      design = factor(
        design,
        levels = c(
          "Grid 700m", "Grid 800m",
          "Cluster (OS)", "Cluster (2 Sig)",
          "Lacework",
          "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both",
          "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both",
          "Two Stage"
        ),
        labels = c(
          "Gr.700", "Gr.800",
          "Cl(OS)", "Cl(2 Sig)",
          "LW",
          "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both", 
          "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both",
          "2 Stage"
        )
      )
    )
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
      ~ stratum + traps,
      labeller = labeller(
        stratum = c("S1" = "Stratum 1",
                    "S2" = "Stratum 2"),
        traps = c("40" = "40 traps",
                  "120" = "120 traps")
      )
    )
    
  } else {
    
    p <- p + facet_wrap(
      ~ stratum,
      labeller = labeller(
        stratum = c("S1" = "Stratum 1",
                    "S2" = "Stratum 2")
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
Combine.plots <- function(p1, p2, p3, global_title = NULL, tag_levels = "A", legend_position = "top", compact = TRUE) {
  
  # ✅ Base combination
  p_comb <- (p1 / p2 / p3) +
    plot_layout(
      guides = "collect",   # ✅ shared legend
      heights = c(1, 1, 1)  # equal heights
    )
  
  # ✅ Compact theme for vertical space
  if (compact) {
    p_comb <- p_comb &
      theme(
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
  }
  
  # ✅ Add labels + title
  p_comb <- p_comb +
    plot_annotation(
      title = global_title,
      tag_levels = tag_levels
    )
  
  return(p_comb)
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

