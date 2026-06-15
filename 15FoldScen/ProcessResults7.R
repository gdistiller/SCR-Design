#June 2026, process sim results for 15 fold scen, 500 reps
#For 40 traps and 120 traps and a large regular area
#filters out reps where se/est > 100
#this version first constructs a big df with all results

library(secr)
library(secrdesign)
library(dplyr)
library(tidyr)
library(ggplot2)

################
#Functions
################

#function to construct the df with all results
construct_df <- function(df, design_name, true_vals) {
  
  # Optional safety check (highly recommended)
  expected_params <- c("D", "L0", "Sig")
  for (p in expected_params) {
    stopifnot(paste0(p, ".x") %in% names(df))
    stopifnot(paste0(p, ".se.x") %in% names(df))
    stopifnot(paste0(p, ".y") %in% names(df))
    stopifnot(paste0(p, ".se.y") %in% names(df))
  }
  
  df %>%
    mutate(sim = row_number()) %>%
    
    # reshape to long
    pivot_longer(
      cols = matches("\\.(x|y)$"),
      names_to = "var",
      values_to = "value"
    ) %>%
    
    # extract stratum
    mutate(
      stratum = sub(".*\\.", "", var),
      stratum = ifelse(stratum == "x", "S1", "S2")
    ) %>%
    
    # remove stratum suffix
    mutate(
      var_nostrat = sub("\\.(x|y)$", "", var)
    ) %>%
    
    # identify type and parameter
    mutate(
      type = ifelse(grepl("\\.se$", var_nostrat), "se", "estimate"),
      param = sub("\\.se$", "", var_nostrat)
    ) %>%
    
    select(sim, param, stratum, type, value) %>%
    
    # reshape back to wide (estimate + se)
    pivot_wider(
      names_from = type,
      values_from = value
    ) %>%
    
    mutate(design = design_name) %>%
    
    # join true values safely
    left_join(true_vals, by = c("param", "stratum"))
}

#function to find rogue values for a single strata
#cannot do both together as diff strata cant be in the same row
#true values provided as D, L0, sigma
#detects NA and infinite values, and for any estimates greater or less than mag fold 
find.rogue <- function(df, mag = 10, true){
  bad.rows <- NULL
  df <- as.data.frame(df)
  colnames(df) <- c("D", "D.se", "L0", "L0.se", "Sig", "Sig.se")
  
  ##identify rows with NA or infinite
  inds1 <- which(is.na(rowSums(df))|is.infinite(rowSums(df)))
  
  ##identify rows where the relative se is v large relative to the est
  RSE.D <- df[,2] / df[,1] 
  RSE.L0 <- df[,4] / df[,3] 
  RSE.Sig <- df[,6] / df[,5]
  
  inds2 <- which(RSE.D > 100 | RSE.L0 > 100 | RSE.Sig > 100)
  
  #identify rows with wild values for any estimates (mag fold)
  inds3 <- which(df$D > true[1]*mag | df$D < true[1]/mag | df$L0 > true[2]*mag | df$L0 < true[2]/mag | df$Sig > true[3]*mag | df$Sig < true[3]/mag)  
  
  bad.rows <- c(inds1, inds2, inds3)
  
  if (length(bad.rows)==0) bad.rows <- NULL
  
  if (!is.null(bad.rows)) {
    bad.rows <- sort(unique(bad.rows)) 
    df <- df[-bad.rows,]
  }
  
  df$row <- row.names(df)
  return("Clean data" = df)
}

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
        grepl("GA4", design) ~ "Optimised 1",
        grepl("GA5", design) ~ "Optimised 2",
        TRUE ~ "Other"
      ),
      
      design_group = factor(
        design_group,
        levels = c("Grid", "Cluster", "Optimised 1", "Optimised 2")
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
          "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both", "GA4 Both (max-min)",
          "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both", "GA5 BothMax"
        ),
        labels = c(
          "Gr.700", "Gr.800",
          "Cl(OS)", "Cl(2 Sig)",
          "LW",
          "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both", "GA4 Both (max-min)",
          "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both", "GA5 BothMax"
        )
      )
    )
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
}

####################################################
#set some global values for plotting and filtering
mag.factor <- 5 
D1.ylim <- 0.25 ; D2.ylim <- 0.015
L1.ylim <- 5 ; L2.ylim <- 0.3
Sig1.ylim <- 500 ; Sig2.ylim <- 10000
se.lim <- 12

true.vals <- matrix(c(
  0.05,  0.05/15,    # D.x, D.y
  2,  2/15,  # L0.x, L0.y
  200, 3000 # Sig.x, Sig.y
), nrow = 3, byrow = TRUE)

rownames(true.vals) <- c("D", "L0", "Sig")
colnames(true.vals) <- c("S1", "S2")

truth_df <- as.data.frame(true.vals) %>%
  tibble::rownames_to_column("param") %>%
  tidyr::pivot_longer(
    cols = -param,
    names_to = "stratum",
    values_to = "truth"
  )
####################################################

#################
#Build data frames
#40 traps first
#################

###########################################################
#Grid designs (incl grids with OS 2G, and clustered grids)
###########################################################

#800 m spacing

load("15FoldScen/Cluster/Sims/40Traps/GridsdResults40.RData")

Grid40.800.red1 <- find.rogue(Grid.800.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid40.800.red2 <- find.rogue(Grid.800.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid40.800.clean <- merge(Grid40.800.red1, Grid40.800.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid800.40 <- construct_df(Grid40.800.clean, design_name = "Grid 800m", true_vals = truth_df)

#700 m spacing (OS)
load("15FoldScen/Cluster/Sims/40Traps/OSGridsResults40b.RData")

Grid40.700.red1 <- find.rogue(Grid.700.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid40.700.red2 <- find.rogue(Grid.700.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid40.700.clean <- merge(Grid40.700.red1, Grid40.700.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid700.40 <- construct_df(Grid40.700.clean, design_name = "Grid 700m", true_vals = truth_df)

#cluster (os spacing)
load("15FoldScen/Cluster/Sims/40Traps/Clusters1b.RData")

Cluster40.os.red1 <- find.rogue(Grid.os.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster40.os.red2 <- find.rogue(Grid.os.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Cluster40.os.clean <- merge(Cluster40.os.red1, Cluster40.os.red2, by = "row", all.x = TRUE, all.y = TRUE)

ClusterOS.40 <- construct_df(Cluster40.os.clean, design_name = "Cluster (OS)", true_vals = truth_df)

#now 2 sigma spacing
load("15FoldScen/Cluster/Sims/40Traps/Clusters2.RData")

Cluster40.2sig.red1 <- find.rogue(Grid.2sig.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster40.2sig.red2 <- find.rogue(Grid.2sig.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Cluster40.2sig.clean <- merge(Cluster40.2sig.red1, Cluster40.2sig.red2, by = "row", all.x = TRUE, all.y = TRUE)

Cluster2Sig.40 <- construct_df(Cluster40.2sig.clean, design_name = "Cluster (2 Sig)", true_vals = truth_df)

################
#GA4 designs
################

load("15FoldScen/Cluster/Sims/40Traps/GA4dResults40.RData")

#S1 pars
GA4.40.S1.red1 <- find.rogue(G4.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.S1.red2 <- find.rogue(G4.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.40.S1.clean <- merge(GA4.40.S1.red1, GA4.40.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.S1.40 <- construct_df(GA4.40.S1.clean, design_name = "GA4 S1", true_vals = truth_df)

#S2 pars
GA4.40.S2.red1 <- find.rogue(G4.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.S2.red2 <- find.rogue(G4.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.40.S2.clean <- merge(GA4.40.S2.red1, GA4.40.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.S2.40 <- construct_df(GA4.40.S2.clean, design_name = "GA4 S2", true_vals = truth_df)

#Avg pars
GA4.40.Avg.red1 <- find.rogue(G4.Avg.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.Avg.red2 <- find.rogue(G4.Avg.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.40.Avg.clean <- merge(GA4.40.Avg.red1, GA4.40.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.Avg.40 <- construct_df(GA4.40.Avg.clean, design_name = "GA4 Avg", true_vals = truth_df)

#Both
GA4.40.Both.red1 <- find.rogue(G4.Both.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.Both.red2 <- find.rogue(G4.Both.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.40.Both.clean <- merge(GA4.40.Both.red1, GA4.40.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.Both.40 <- construct_df(GA4.40.Both.clean, design_name = "GA4 Both", true_vals = truth_df)

#Both (max-min)
GA4.40.BothMax.red1 <- find.rogue(G4.BothMax.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.BothMax.red2 <- find.rogue(G4.BothMax.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.40.BothMax.clean <- merge(GA4.40.BothMax.red1, GA4.40.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.BothMax.40 <- construct_df(GA4.40.BothMax.clean, design_name = "GA4 Both (max-min)", true_vals = truth_df)

#############
#GA5 designs
#############

load("15FoldScen/Cluster/Sims/40Traps/GA5dResults40.RData")

#S1 pars
GA5.40.S1.red1 <- find.rogue(G5.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.S1.red2 <- find.rogue(G5.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.40.S1.clean <- merge(GA5.40.S1.red1, GA5.40.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.S1.40 <- construct_df(GA5.40.S1.clean, design_name = "GA5 S1", true_vals = truth_df)

#S2 pars
GA5.40.S2.red1 <- find.rogue(G5.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.S2.red2 <- find.rogue(G5.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.40.S2.clean <- merge(GA5.40.S2.red1, GA5.40.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.S2.40 <- construct_df(GA5.40.S2.clean, design_name = "GA5 S2", true_vals = truth_df)

#Avg pars
GA5.40.Avg.red1 <- find.rogue(G5.Avg.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.Avg.red2 <- find.rogue(G5.Avg.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.40.Avg.clean <- merge(GA5.40.Avg.red1, GA5.40.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.Avg.40 <- construct_df(GA5.40.Avg.clean, design_name = "GA5 Avg", true_vals = truth_df)

#Both
GA5.40.Both.red1 <- find.rogue(G5.Both.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.Both.red2 <- find.rogue(G5.Both.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.40.Both.clean <- merge(GA5.40.Both.red1, GA5.40.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.Both.40 <- construct_df(GA5.40.Both.clean, design_name = "GA5 Both", true_vals = truth_df)

#Both (max-min)
GA5.40.BothMax.red1 <- find.rogue(G5.BothMax.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.BothMax.red2 <- find.rogue(G5.BothMax.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.40.BothMax.clean <- merge(GA5.40.BothMax.red1, GA5.40.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.BothMax.40 <- construct_df(GA5.40.BothMax.clean, design_name = "GA5 BothMax", true_vals = truth_df)

#########
#Lacework
#########

load("15FoldScen/Cluster/Sims/40Traps/LWResults40.RData")

LW40.red1 <- find.rogue(LW.40.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LW40.red2 <- find.rogue(LW.40.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
LW40.clean <- merge(LW40.red1, LW40.red2, by = "row", all.x = TRUE, all.y = TRUE)

LW.40 <- construct_df(LW40.clean, design_name = "Lacework", true_vals = truth_df)


all_results_40 <- bind_rows(
  Grid800.40,
  Grid700.40,
  ClusterOS.40,
  Cluster2Sig.40,
  GA4.S1.40,
  GA4.S2.40,
  GA4.Avg.40,
  GA4.Both.40,
  GA4.BothMax.40,
  GA5.S1.40,
  GA5.S2.40,
  GA5.Avg.40,
  GA5.Both.40,
  GA5.BothMax.40,
  LW.40)

#####################################

#################
#120 traps first
#################

###########################################################
#Grid designs (incl grids with OS 2G, and clustered grids)
###########################################################

#800 m spacing

load("15FoldScen/Cluster/Sims/120Traps/Grids800Results120.RData")

Grid120.800.red1 <- find.rogue(Grid.800.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid120.800.red2 <- find.rogue(Grid.800.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid120.800.clean <- merge(Grid120.800.red1, Grid120.800.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid800.120 <- construct_df(Grid120.800.clean, design_name = "Grid 800m", true_vals = truth_df)

#700 m spacing (OS)
load("15FoldScen/Cluster/Sims/120Traps/Grids700Results120.RData")

Grid120.700.red1 <- find.rogue(Grid.700.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid120.700.red2 <- find.rogue(Grid.700.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid120.700.clean <- merge(Grid120.700.red1, Grid120.700.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid700.120 <- construct_df(Grid120.700.clean, design_name = "Grid 700m", true_vals = truth_df)

#cluster (os spacing)
load("15FoldScen/Cluster/Sims/120Traps/Clusters120a.RData")

Cluster120.os.red1 <- find.rogue(Grid.os.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster120.os.red2 <- find.rogue(Grid.os.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Cluster120.os.clean <- merge(Cluster120.os.red1, Cluster120.os.red2, by = "row", all.x = TRUE, all.y = TRUE)

ClusterOS.120 <- construct_df(Cluster120.os.clean, design_name = "Cluster (OS)", true_vals = truth_df)

#now 2 sigma spacing
load("15FoldScen/Cluster/Sims/120Traps/Clusters120b.RData")

Cluster120.2sig.red1 <- find.rogue(Grid.2sig.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster120.2sig.red2 <- find.rogue(Grid.2sig.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Cluster120.2sig.clean <- merge(Cluster120.2sig.red1, Cluster120.2sig.red2, by = "row", all.x = TRUE, all.y = TRUE)

Cluster2Sig.120 <- construct_df(Cluster120.2sig.clean, design_name = "Cluster (2 Sig)", true_vals = truth_df)

################
#GA4 designs
################

#S1 pars
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120a.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120a2.RData")
G4.S1.results <- rbind(G4.S1.resultsa, G4.S1.resultsa2)

GA4.120.S1.red1 <- find.rogue(G4.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.S1.red2 <- find.rogue(G4.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.120.S1.clean <- merge(GA4.120.S1.red1, GA4.120.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.S1.120 <- construct_df(GA4.120.S1.clean, design_name = "GA4 S1", true_vals = truth_df)

#S2 pars
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120b.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120b2.RData")
G4.S2.results <- rbind(G4.S2.resultsb, G4.S2.resultsb2)

GA4.120.S2.red1 <- find.rogue(G4.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.S2.red2 <- find.rogue(G4.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.120.S2.clean <- merge(GA4.120.S2.red1, GA4.120.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.S2.120 <- construct_df(GA4.120.S2.clean, design_name = "GA4 S2", true_vals = truth_df)

#Avg pars
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120C.RData")

GA4.120.Avg.red1 <- find.rogue(G4.Avg.resultsc[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.Avg.red2 <- find.rogue(G4.Avg.resultsc[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.120.Avg.clean <- merge(GA4.120.Avg.red1, GA4.120.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.Avg.120 <- construct_df(GA4.120.Avg.clean, design_name = "GA4 Avg", true_vals = truth_df)

#Both
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120d.RData")

GA4.120.Both.red1 <- find.rogue(G4.Both.resultsd[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.Both.red2 <- find.rogue(G4.Both.resultsd[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.120.Both.clean <- merge(GA4.120.Both.red1, GA4.120.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.Both.120 <- construct_df(GA4.120.Both.clean, design_name = "GA4 Both", true_vals = truth_df)

#Both (max-min)
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120e.RData")
GA4.120.BothMax.red1 <- find.rogue(G4.BothMax.resultse[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.BothMax.red2 <- find.rogue(G4.BothMax.resultse[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA4.120.BothMax.clean <- merge(GA4.120.BothMax.red1, GA4.120.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA4.BothMax.120 <- construct_df(GA4.120.BothMax.clean, design_name = "GA4 Both (max-min)", true_vals = truth_df)

#############
#GA5 designs
#############

#S1 pars
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120a.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120a2.RData")
G5.S1.results <- rbind(G5.S1.resultsa, G5.S1.resultsa2)

GA5.120.S1.red1 <- find.rogue(G5.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.S1.red2 <- find.rogue(G5.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.120.S1.clean <- merge(GA5.120.S1.red1, GA5.120.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.S1.120 <- construct_df(GA5.120.S1.clean, design_name = "GA5 S1", true_vals = truth_df)

#S2 pars
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120b.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120b2.RData")
G5.S2.results <- rbind(G5.S2.resultsb, G5.S2.resultsb2)

GA5.120.S2.red1 <- find.rogue(G5.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.S2.red2 <- find.rogue(G5.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.120.S2.clean <- merge(GA5.120.S2.red1, GA5.120.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.S2.120 <- construct_df(GA5.120.S2.clean, design_name = "GA5 S2", true_vals = truth_df)

#Avg pars
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120c.RData")

GA5.120.Avg.red1 <- find.rogue(G5.Avg.resultsc[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.Avg.red2 <- find.rogue(G5.Avg.resultsc[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.120.Avg.clean <- merge(GA5.120.Avg.red1, GA5.120.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.Avg.120 <- construct_df(GA5.120.Avg.clean, design_name = "GA5 Avg", true_vals = truth_df)

#Both
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120d.RData")

GA5.120.Both.red1 <- find.rogue(G5.Both.resultsd[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.Both.red2 <- find.rogue(G5.Both.resultsd[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.120.Both.clean <- merge(GA5.120.Both.red1, GA5.120.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.Both.120 <- construct_df(GA5.120.Both.clean, design_name = "GA5 Both", true_vals = truth_df)

#Both (max-min)
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120e.RData")

GA5.120.BothMax.red1 <- find.rogue(G5.BothMax.resultse[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.BothMax.red2 <- find.rogue(G5.BothMax.resultse[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
GA5.120.BothMax.clean <- merge(GA5.120.BothMax.red1, GA5.120.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)

GA5.BothMax.120 <- construct_df(GA5.120.BothMax.clean, design_name = "GA5 BothMax", true_vals = truth_df)

#########
#Lacework
#########

load("15FoldScen/Cluster/Sims/120Traps/LWResults120.RData")

LW120.red1 <- find.rogue(LW.120.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LW120.red2 <- find.rogue(LW.120.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
LW120.clean <- merge(LW120.red1, LW120.red2, by = "row", all.x = TRUE, all.y = TRUE)

LW.120 <- construct_df(LW120.clean, design_name = "Lacework", true_vals = truth_df)

################################################################

all_results_120 <- bind_rows(
  Grid800.120,
  Grid700.120,
  ClusterOS.120,
  Cluster2Sig.120,
  GA4.S1.120,
  GA4.S2.120,
  GA4.Avg.120,
  GA4.Both.120,
  GA4.BothMax.120,
  GA5.S1.120,
  GA5.S2.120,
  GA5.Avg.120,
  GA5.Both.120,
  GA5.BothMax.120,
  LW.120)

Allresults <- list("40 traps" = all_results_40, "120 traps" = all_results_120)
save(Allresults, file = "15FoldScen/AllResults.RData")

########################################################################
#Plots
library(ggplot2)
library(dplyr)
library(tidyr)

##############################################
#load results

load("15FoldScen/AllResults.RData")

AllResults40 <- Allresults[[1]]
AllResults120 <- Allresults[[2]]

###############################################
#create summaries
summary.40 <- make.summary(AllResults40)
summary.120 <- make.summary(AllResults120)

summary.40$traps  <- "40"
summary.120$traps <- "120"
summary.combined <- dplyr::bind_rows(summary.40, summary.120)

#RB plot, 40 traps
RB.plot(summary.40, "Relative Bias (40 traps)", ylims = c(0, 0.5))

#RB plot, 120 traps
RB.plot(summary.120, "Relative Bias (120 traps)", ylims = c(0, 0.5))

#RB.plot 2 operates on a combined df to nest trap lvl within stratum
#x axis labels are too squashed
RB.plot2(summary.combined, "Relative Bias (40 traps)", ylims = c(0, 0.5))

#Precision plot, 40 traps
Prec.plot(summary.40, "Relative SE (40 traps)", ylims = c(0, 1.5))

#Precision plot, 120 traps
Prec.plot(summary.120, "Relative SE (120 traps)", ylims = c(0, 1.5))

#Coverage, 40 traps
Coverage.plot(summary.40, "Coverage of 95% Confidence Intervals (40 traps)", ylims = c(0.8,1))

#Coverage, 120 traps
Coverage.plot(summary.120, "Coverage of 95% Confidence Intervals (120 traps)", ylims = c(0.8,1))


