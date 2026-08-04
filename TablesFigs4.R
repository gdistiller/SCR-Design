#Aug 2026, new plots (scatter) 
#uses AllResults4.RData

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggrepel)

#################################################################
#figs created in ggplot and saved
#needs make.summary, and new scatter plot fn

load("15FoldScen/AllResults4.RData")
AllResults40 <- Allresults_4[[1]]
AllResults120 <- Allresults_4[[2]]

#control order of design factor
AllResults40 <- AllResults40 %>%
  mutate(
    design = factor(
      design,
      levels = c(
        "Grid (OS)", "Grid 800m",
        "Cluster (OS)", "Cluster (2 Sig)",
        "Lacework", "Lacework (F)",
        "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both",
        "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both",
        "Two Stage"
      )
    )
  )

AllResults120 <- AllResults120 %>%
  mutate(
    design = factor(
      design,
      levels = c(
        "Grid (OS)", "Grid 800m",
        "Cluster (OS)", "Cluster (2 Sig)",
        "Lacework", "Lacework (F)",
        "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both",
        "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both",
        "Two Stage"
      )
    )
  )


#create summaries
summary.40 <- make.summary(AllResults40)
summary.120 <- make.summary(AllResults120)
summary.40$traps  <- "40"
summary.120$traps <- "120"
summary.combined <- dplyr::bind_rows(summary.40, summary.120)
summary.combined$traps <- factor(summary.combined$traps, levels = c("40", "120"))

#Define scale
# ----------------------------------------------------
# Common colour palette and legend
# ----------------------------------------------------

design_cols <- c(
  "Grid"     = "#0072B2",
  "Cluster"  = "#D55E00",
  "Lacework" = "#009E73",
  "GA4"      = "#CC79A7",
  "GA5"      = "#E69F00",
  "2 Stage"  = "#000000"
)

design_colour_scale <- scale_colour_manual(
  values = design_cols,
  name = "Design family",
  drop = FALSE
)

design_guide <- guides(
  colour = guide_legend(
    nrow = 1,
    byrow = TRUE,
    override.aes = list(
      shape = 16,
      size = 3,
      linewidth = 0,
      linetype = 0
    )
  )
)

#short labels for the plots
design_key <- tibble::tribble(
  ~design,       ~short_label,
  "Grid (OS)",   "Gr1",
  "Grid 800",    "Gr2",
  
  "Cl (OS)",     "Cl1",
  "Cl (2 Sig)",  "Cl2",
  
  "LW",          "LW1",
  "LW (f)",      "LW2",
  
  "GA4 G1",      "GA4-1",
  "GA4 G2",      "GA4-2",
  "GA4 Avg",     "GA4-3",
  "GA4 Both",    "GA4-4",
  
  "GA5 G1",      "GA5-1",
  "GA5 G2",      "GA5-2",
  "GA5 Avg",     "GA5-3",
  "GA5 Both",    "GA5-4",
  
  "2 Stage",     "TS"
)

#new scatter plots for RB and RSE
D.RB40 <- make_metric_scatter(summary.combined, metric = "RB", trap_level = "40",
                    parameter = "Density")

D.RB120 <- make_metric_scatter(summary.combined, metric = "RB", trap_level = "120",
                              parameter = "Density")

D.RSE40 <- make_metric_scatter(summary.combined, metric = "RSE_emp", trap_level = "40",
                              parameter = "Density")

D.RSE120 <- make_metric_scatter(summary.combined, metric = "RSE_emp", trap_level = "120",
                               parameter = "Density")

#new coverage plot
D.Cov40 <- Coverage.plot(summary.combined,
                     plot.title = "40 traps",
                     trap_level = "40",
                     param_select = "Density",
                     ylims=c(0.8,1))

D.Cov120 <- Coverage.plot(summary.combined,
                         plot.title = "120 traps",
                         trap_level = "120",
                         param_select = "Density",
                         ylims=c(0.8,1))


#combine
plot.D2 <- Combine.performance.plots(
  D.RB40,
  D.RB120,
  D.RSE40,
  D.RSE120,
  D.Cov40,
  D.Cov120,
  global_title = "Performance for Density"
)

ggsave("15FoldScen/figures/D2.pdf", plot = plot.D2, width = 7, height = 9)

#sigma
#new scatter plots for RB and RSE
Sig.RB40 <- make_metric_scatter(summary.combined, metric = "RB", trap_level = "40",
                              parameter = "sigma")

Sig.RB120 <- make_metric_scatter(summary.combined, metric = "RB", trap_level = "120",
                               parameter = "sigma")

Sig.RSE40 <- make_metric_scatter(summary.combined, metric = "RSE_emp", trap_level = "40",
                               parameter = "sigma")

Sig.RSE120 <- make_metric_scatter(summary.combined, metric = "RSE_emp", trap_level = "120",
                                parameter = "sigma")

#coverage plots
Sig.Cov40 <- Coverage.plot(summary.combined,
                         plot.title = "40 traps",
                         trap_level = "40",
                         param_select = "sigma",
                         ylims=c(0.8,1))

Sig.Cov120 <- Coverage.plot(summary.combined,
                          plot.title = "120 traps",
                          trap_level = "120",
                          param_select = "sigma",
                          ylims=c(0.8,1))


#combine
plot.Sig2 <- Combine.performance.plots(
  Sig.RB40,
  Sig.RB120,
  Sig.RSE40,
  Sig.RSE120,
  Sig.Cov40,
  Sig.Cov120,
  global_title = expression("Performance of " * sigma)
)

ggsave("15FoldScen/figures/Sig2.pdf", plot = plot.Sig2, width = 7, height = 9)

#L0
#new scatter plots for RB and RSE
L0.RB40 <- make_metric_scatter(summary.combined, metric = "RB", trap_level = "40",
                                parameter = "lambda[0]")

L0.RB120 <- make_metric_scatter(summary.combined, metric = "RB", trap_level = "120",
                                 parameter = "lambda[0]")

L0.RSE40 <- make_metric_scatter(summary.combined, metric = "RSE_emp", trap_level = "40",
                                 parameter = "lambda[0]")

L0.RSE120 <- make_metric_scatter(summary.combined, metric = "RSE_emp", trap_level = "120",
                                  parameter = "lambda[0]")

#coverage plots
L0.Cov40 <- Coverage.plot(summary.combined,
                           plot.title = "40 traps",
                           trap_level = "40",
                           param_select = "lambda[0]",
                           ylims=c(0.8,1))

L0.Cov120 <- Coverage.plot(summary.combined,
                            plot.title = "120 traps",
                            trap_level = "120",
                            param_select = "lambda[0]",
                            ylims=c(0.8,1))


#combine
plot.L02 <- Combine.performance.plots(
  L0.RB40,
  L0.RB120,
  L0.RSE40,
  L0.RSE120,
  L0.Cov40,
  L0.Cov120,
  global_title = expression("Performance of " * lambda[0])
)

ggsave("15FoldScen/figures/L02.pdf", plot = plot.L02, width = 7, height = 9)

