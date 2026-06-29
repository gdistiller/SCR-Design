#For tables, this script uses kableextra to construct tables in R
#then the latex code behind the table is saved to a .tex file
#lastly the code is included in oOverleaf using: \input{tables/table1.tex}

library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# table of parameter estimates per strata

true.vals <- matrix(c(
  0.05,  0.05/15,    # D.x, D.y
  2,  2/15,  # L0.x, L0.y
  200, 3000 # Sig.x, Sig.y
), nrow = 3, byrow = TRUE)


rownames(true.vals) <- c("$D$", "$\\lambda_0$", "$\\sigma$")
colnames(true.vals) <- c("S1", "S2")

#format the values before passing to kableextra

true.vals_fmt <- true.vals

# format each row differently
true.vals_fmt[1, ] <- formatC(true.vals[1, ], format = "f", digits = 3)        # integer
true.vals_fmt[2, ] <- formatC(true.vals[2, ], format = "f", digits = 3)
true.vals_fmt[3, ] <- formatC(true.vals[3, ], format = "d")

tab1 <- kable(true.vals_fmt,
        format = "latex",
        booktabs = TRUE,
        escape = FALSE)


writeLines(tab1, "15FoldScen/tables/tab1.tex")


#################################################################
#figs created in ggplot and saved
#needs make.summary, Metric.plot and Combine.plots fns

#load results and set things up
load("15FoldScen/AllResults3.RData")
AllResults40 <- Allresults_3[[1]]
AllResults120 <- Allresults_3[[2]]
AllResults40$traps  <- "40"
AllResults120$traps <- "120"
AllResults.combined <- bind_rows(AllResults40, AllResults120)
AllResults.combined$traps <- factor(AllResults.combined$traps, levels = c("40", "120"))
AllResults.combined$relSE <- AllResults.combined$se / AllResults.combined$truth

#control order of design factor
AllResults.combined <- AllResults.combined %>%
  mutate(
    design = factor(
      design,
      levels = c(
        "Grid 700m", "Grid 800m",
        "Cluster (OS)", "Cluster (2 Sig)",
        "Lacework",
        "GA4 S1", "GA4 S2", "GA4 Avg", "GA4 Both",
        "GA5 S1", "GA5 S2", "GA5 Avg", "GA5 Both",
        "Two Stage"
      )
    )
  )

design_groups <- list(
  "Grids with 800 m spacing" = c("Grid 800m"),
  "Grids with 700 m spacing" = c("Grid 700m"),
  "Clusters - OS spacing" = c("Cluster (OS)"),
  "Clusters - 2 sigma spacing" = c("Cluster (2 Sig)"),
  "GA4 - S1 values" = c("GA4 S1"),
  "GA4 - S2 values" = c("GA4 S2"),
  "GA4 - Avg values" = c("GA4 Avg"),
  "GA4 - Both values" = c("GA4 Both"),
  "GA5 - S1 values" = c("GA5 S1"),
  "GA5 - S2 values" = c("GA5 S2"),
  "GA5 - Avg values" = c("GA5 Avg"),
  "GA5 - Both values" = c("GA5 Both"),
  "Lacework" = c("Lacework"),
  "Two Stage" = c("2 Stage")
  )

params <- c("D", "L0", "Sig")

#create summaries
summary.40 <- make.summary(AllResults40)
summary.120 <- make.summary(AllResults120)
summary.40$traps  <- "40"
summary.120$traps <- "120"
summary.combined <- dplyr::bind_rows(summary.40, summary.120)
summary.combined$traps <- factor(summary.combined$traps, levels = c("40", "120"))

#plots for D, 1st with trap lvls plotted together
D.RB <- Metric.plot(summary.combined,
                    plot.title = "RB (Density)",
                    param_select = "Density",
                    metric = "RB",
                    facet_traps = FALSE,
                    ylims=c(-0.1,1))

#Rel SE plots (model-based)
D.RSE <- Metric.plot(summary.combined,
                     plot.title = "RSE (Density)",
                     param_select = "Density",
                     metric = "RSE",
                     facet_traps = FALSE,
                     ylims=c(0,1))

#coverage plots
D.cov <- Metric.plot(summary.combined,
                     plot.title = "Coverage (Density)",
                     param_select = "Density",
                     metric = "coverage",
                     facet_traps = FALSE,
                     ylims=c(0.75,1))

#again with traps faceted
D.RB2 <- Metric.plot(summary.combined,
                     plot.title = NULL,
                     param_select = "Density",
                     metric = "RB",
                     facet_traps = TRUE,
                     ylims=c(-0.1,1))

#Rel SE plots (model-based)
D.RSE2 <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "Density",
                      metric = "RSE",
                      facet_traps = TRUE,
                      ylims=c(0,1))

#coverage plots
D.cov2 <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "Density",
                      metric = "coverage",
                      facet_traps = TRUE,
                      ylims=c(0.75,1))

#combine in two different ways
plot.D1 <- Combine.plots(D.RB, D.RSE, D.cov,
              global_title = "Performance for Density"
)

plot.D2 <-Combine.plots(D.RB2, D.RSE2, D.cov2,
              global_title = "Performance for Density"
)

ggsave("15FoldScen/figures/D1.pdf", plot = plot.D1, width = 7, height = 9)
ggsave("15FoldScen/figures/D2.pdf", plot = plot.D2, width = 7, height = 9)

#plots for L0
L0.RB <- Metric.plot(summary.combined,
                    plot.title = NULL,
                    param_select = "lambda[0]",
                    metric = "RB",
                    facet_traps = TRUE,
                    ylims=c(-0.1,1))

#Rel SE plots (model-based)
L0.RSE <- Metric.plot(summary.combined,
                     plot.title = NULL,
                     param_select = "lambda[0]",
                     metric = "RSE",
                     facet_traps = TRUE,
                     ylims=c(0,1))

#coverage plots
L0.cov <- Metric.plot(summary.combined,
                     plot.title = NULL,
                     param_select = "lambda[0]",
                     metric = "coverage",
                     facet_traps = TRUE,
                     ylims=c(0.8,1))

#combine
plot.L01 <- Combine.plots(L0.RB, L0.RSE, L0.cov,
                         global_title = expression("Performance of "~lambda[0])
)
ggsave("15FoldScen/figures/L01.pdf", plot = plot.L01, width = 7, height = 9)

#plots for sigma, trap lvls plotted together
sig.RB <- Metric.plot(summary.combined,
                     plot.title = NULL,
                     param_select = "sigma",
                     metric = "RB",
                     facet_traps = TRUE,
                     ylims=c(-0.1,0.5))

#Rel SE plots (model-based)
sig.RSE <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "sigma",
                      metric = "RSE",
                      facet_traps = TRUE,
                      ylims=c(0,0.5))

#coverage plots
sig.cov <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "sigma",
                      metric = "coverage",
                      facet_traps = TRUE,
                      ylims=c(0.8,1))

#combine in two different ways
plot.sig1 <- Combine.plots(sig.RB, sig.RSE, sig.cov,
                          global_title = expression("Performance of "~ sigma)
)
ggsave("15FoldScen/figures/Sig1.pdf", plot = plot.sig1, width = 7, height = 9)


###################################
#a plot of different realisations, for 40 traps
#objs with various proposed designs pasted into dir below
setwd("~/Git/SCR-Design/15FoldScen/Cluster/ProposedDesigns")
setwd("~/Documents/Git/SCR-Design/15FoldScen/Cluster/ProposedDesigns")

load("SCRObjs.RData")

load("GridDesigns40.RData")
load("LWdesigns.RData")
load("GADesigns40.RData")
load("GA2StageDesigns.RData") 

#extract coordinates and put all in one df

# extract designs
sys1 <- as.data.frame(grid.designs.40$`800 m`[[2]]) %>%
  mutate(design = "Grid 800")

sys2 <- as.data.frame(grid.designs.40$`Cluster (opt)`[[1]]) %>%
  mutate(design = "Cluster (opt)")

sys3 <- as.data.frame(grid.designs.40$`Cluster (2 sig)`[[1]]) %>%
  mutate(design = "Cluster (2 sig)")

sys4 <- as.data.frame(lwlist$`40 traps`) %>%
  mutate(design = "Lacework")

ga41 <- as.data.frame(GA.designs.40$G4S1[[1]]) %>%
  mutate(design = "GA4 S1")

ga42 <- as.data.frame(GA.designs.40$G4S2[[1]]) %>%
  mutate(design = "GA4 S2")

ga43 <- as.data.frame(GA.designs.40$G4Avg[[1]]) %>%
  mutate(design = "GA4 Avg")

ga44 <- as.data.frame(GA.designs.40$G4Both[[1]]) %>%
  mutate(design = "GA4 Both")

ga51 <- as.data.frame(GA.designs.40$G5S1[[1]]) %>%
  mutate(design = "GA5 S1")

ga52 <- as.data.frame(GA.designs.40$G5S2[[1]]) %>%
  mutate(design = "GA5 S2")

ga53 <- as.data.frame(GA.designs.40$G5Avg[[1]]) %>%
  mutate(design = "GA5 Avg")

ga54 <- as.data.frame(GA.designs.40$G5Both[[1]]) %>%
  mutate(design = "GA5 Both")

ga2 <- as.data.frame(GA2StageDesigns$`40 traps`[[1]]) %>%
  mutate(design = "Two Stage")

#first deal with systematic designs
sys.40 <- bind_rows(sys1, sys2, sys3, sys4)

#set labels to parse properly
sys.40 <- sys.40 %>%
  mutate(
    design_label = case_when(
      design == "Grid 800" ~ "Grid~800",
      design == "Cluster (opt)" ~ "Cluster~'(OS)'",
      design == "Cluster (2 sig)" ~ "Cluster~(2*sigma)",
      design == "Lacework" ~ "Lacework"
    ),
    design_label = factor(design_label, levels = c(
      "Grid~800",
      "Cluster~'(OS)'",
      "Cluster~(2*sigma)",
      "Lacework"
      ))
  )

msk <- SCR.objs$Full[[1]]
msk.red <- SCR.objs$Full[[2]]
#msk.red used to  show realisations that excl the buffer

#use fn
sys40.plot <- plot.design(sys.40, msk, view = "full", ndim2 = 4, title = "Systematic trap configurations (40 traps)")

setwd("~/Git/SCR-Design")
ggsave("15FoldScen/figures/sys40.pdf", plot = sys40.plot, width = 9, height = 6)

################################################################################
#now optimised designs plotted on full extent
#first just GA4 or GA5
GA4.40 <- bind_rows(ga41, ga42, ga43, ga44)

GA4.40 <- GA4.40 %>%
  mutate(
    design_label = case_when(
      design == "GA4 S1" ~ "GA[4]~S[1]",
      design == "GA4 S2" ~ "GA[4]~S[2]",
      design == "GA4 Avg" ~ "GA[4]~Avg",
      design == "GA4 Both" ~ "GA[4]~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[4]~S[1]",
      "GA[4]~S[2]",
      "GA[4]~Avg",
      "GA[4]~Both")
      )
  )

GA4.40.plot <- plot.design(GA4.40, msk,  view = "full", ndim2 = 4, title.expr = expression(GA[4]~"trap configurations (40 traps)"))

#GA5
GA5.40 <- bind_rows(ga51, ga52, ga53, ga54)

GA5.40 <- GA5.40 %>%
  mutate(
    design_label = case_when(
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[5]~S[1]",
      "GA[5]~S[2]",
      "GA[5]~Avg",
      "GA[5]~Both"
    ))
  )

GA5.40.plot <- plot.design(GA5.40, msk, view = "full", ndim2 = 4, title.expr = expression(GA[5]~"trap configurations (40 traps)"))

ggsave("15FoldScen/figures/GA440.pdf", plot = GA4.40.plot, width = 9, height = 6)
ggsave("15FoldScen/figures/GA540.pdf", plot = GA5.40.plot, width = 9, height = 6)

#and both sets in one plot
#current fn doesn't support this, needs adjustment
opt.40 <- bind_rows(ga41, ga42, ga43, ga44,
                    ga51, ga52, ga53, ga54)

opt.40 <- opt.40 %>%
  mutate(
    design_label = case_when(
      design == "GA4 S1" ~ "GA[4]~S[1]",
      design == "GA4 S2" ~ "GA[4]~S[2]",
      design == "GA4 Avg" ~ "GA[4]~Avg",
      design == "GA4 Both" ~ "GA[4]~Both",
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[4]~S[1]",
      "GA[5]~S[1]",
      "GA[4]~S[2]",
      "GA[5]~S[2]",
      "GA[4]~Avg",
      "GA[5]~Avg",
      "GA[4]~Both",
      "GA[5]~Both"
    ))
  )

opt.40.plot <- plot.design(opt.40, msk, view = "full", ndim1 = 2, ndim2 = 4, title.expr = expression(GA[4]~"and"~GA[5]~"trap configurations (40 traps)"))

ggsave("15FoldScen/figures/Opt40.pdf", plot = opt.40.plot, width = 270, height = 190, units = "mm")

#and both sets in one plot
opt2.40 <- bind_rows(ga41, ga42, ga43, ga44,
                    ga51, ga52, ga53, ga54,
                    ga2)

opt2.40 <- opt2.40 %>%
  mutate(
    design_label = case_when(
      design == "GA4 S1" ~ "GA[4]~S[1]",
      design == "GA4 S2" ~ "GA[4]~S[2]",
      design == "GA4 Avg" ~ "GA[4]~Avg",
      design == "GA4 Both" ~ "GA[4]~Both",
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[4]~S[1]",
      "GA[5]~S[1]",
      "GA[4]~S[2]",
      "GA[5]~S[2]",
      "GA[4]~Avg",
      "GA[5]~Avg",
      "GA[4]~Both",
      "GA[5]~Both",
      "Two~Stage"
    ))
  )

opt2.40.plot <- plot.design(opt2.40, msk, view = "full", ndim1 = 3, ndim2 = 4, title.expr = "Optimised trap configurations (40 traps)")

ggsave("15FoldScen/figures/Opt240.pdf", plot = opt2.40.plot, width = 7, height = 9.5)

#all together, build each row by itself and then combine

#put GA5 and 2 stage together
GA2.40 <- bind_rows(ga51, ga52, ga53, ga54,
                     ga2)

GA2.40 <- GA2.40 %>%
  mutate(
    design_label = case_when(
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[5]~S[1]",
      "GA[5]~S[2]",
      "GA[5]~Avg",
      "GA[5]~Both",
      "Two~Stage"
    ))
  )

#now adding dummy panels to get proper alignment
sys.pad.40 <- bind_rows(
  sys.40,
  tibble(
    x = 0, y = 0,
    design = "dummy",
    design_label = "dummy"
  )
)

levels_all_sys <- c(levels(sys.40$design_label), "dummy")
sys.pad.40$design_label <- factor(sys.pad.40$design_label, levels = levels_all_sys)

ga4.pad.40 <- bind_rows(
  GA4.40,
  tibble(
    x = 0, y = 0,
    design = "dummy",
    design_label = "dummy"
  )
)

levels_all_ga4 <- c(levels(GA4.40$design_label), "dummy")
ga4.pad.40$design_label <- factor(ga4.pad.40$design_label, levels = levels_all_ga4)

p.GA2.40 <- plot.design(GA2.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, 
                        point.size = 0.5, title.expr = NULL)
p.syspad.40 <- plot.design(sys.pad.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, title.expr = NULL, 
                           point.size = 0.5, buffer.prop = 0.05, levels_all = levels_all_sys)
p.GA4pad.40 <- plot.design(ga4.pad.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, 
                           point.size = 0.5, title.expr = NULL, levels_all = levels_all_ga4)

all.40.plot <- Combine.layoutplots(p.syspad.40, p.GA4pad.40,p.GA2.40)

setwd("~/Git/SCR-Design")
ggsave("15FoldScen/figures/All40.pdf", plot = all.40.plot, width = 270, height = 190, units = "mm")

###############
#for 120 traps
#version xxb is from ngen = 500
###############
setwd("~/Git/SCR-Design/15FoldScen/Cluster/ProposedDesigns")
load("GridDesigns120.RData")
load("GADesigns120.RData")
load("LWdesigns.RData")
load("GA2StageDesignsb.RData") 

#extract coordinates and put all in one df

# extract designs
sys1b <- as.data.frame(grid.designs.120$`800 m`[[2]]) %>%
  mutate(design = "Grid 800")

sys2b <- as.data.frame(grid.designs.120$`Cluster (opt)`[[1]]) %>%
  mutate(design = "Cluster (opt)")

sys3b <- as.data.frame(grid.designs.120$`Cluster (2 sig)`[[1]]) %>%
  mutate(design = "Cluster (2 sig)")

sys4b <- as.data.frame(lwlist$`120 traps`) %>%
  mutate(design = "Lacework")

ga41b <- as.data.frame(GA.designs.120$G4S1[[1]]) %>%
  mutate(design = "GA4 S1")

ga42b <- as.data.frame(GA.designs.120$G4S2[[1]]) %>%
  mutate(design = "GA4 S2")

ga43b <- as.data.frame(GA.designs.120$G4Avg[[1]]) %>%
  mutate(design = "GA4 Avg")

ga44b <- as.data.frame(GA.designs.120$G4Both[[1]]) %>%
  mutate(design = "GA4 Both")

ga51b <- as.data.frame(GA.designs.120$G5S1[[1]]) %>%
  mutate(design = "GA5 S1")

ga52b <- as.data.frame(GA.designs.120$G5S2[[1]]) %>%
  mutate(design = "GA5 S2")

ga53b <- as.data.frame(GA.designs.120$G5Avg[[1]]) %>%
  mutate(design = "GA5 Avg")

ga54b <- as.data.frame(GA.designs.120$G5Both[[1]]) %>%
  mutate(design = "GA5 Both")

ga2b <- as.data.frame(GA2StageDesigns$`120 traps`[[1]]) %>%
  mutate(design = "Two Stage")

#first deal with systematic designs
sys.120 <- bind_rows(sys1b, sys2b, sys3b, sys4b)

#set labels to parse properly
sys.120 <- sys.120 %>%
  mutate(
    design_label = case_when(
      design == "Grid 800" ~ "Grid~800",
      design == "Cluster (opt)" ~ "Cluster~'(OS)'",
      design == "Cluster (2 sig)" ~ "Cluster~(2*sigma)",
      design == "Lacework" ~ "Lacework"
    ),
    design_label = factor(design_label, levels = c(
      "Grid~800",
      "Cluster~'(OS)'",
      "Cluster~(2*sigma)",
      "Lacework"
    ))
  )

msk <- SCR.objs$Full[[1]]

#use fn
sys120.plot <- plot.design(sys.120, msk, view = "full", title = "Systematic trap configurations (120 traps)")

setwd("~/Git/SCR-Design")
ggsave("15FoldScen/figures/sys120.pdf", plot = sys120.plot, width = 9, height = 6)

#######################################################################################
#now optimised designs plotted on full extent
#first just GA4 or GA5
GA4.120 <- bind_rows(ga41b, ga42b, ga43b, ga44b)

GA4.120 <- GA4.120 %>%
  mutate(
    design_label = case_when(
      design == "GA4 S1" ~ "GA[4]~S[1]",
      design == "GA4 S2" ~ "GA[4]~S[2]",
      design == "GA4 Avg" ~ "GA[4]~Avg",
      design == "GA4 Both" ~ "GA[4]~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[4]~S[1]",
      "GA[4]~S[2]",
      "GA[4]~Avg",
      "GA[4]~Both")
    )
  )

GA4.120.plot <- plot.design(GA4.120, msk, view = "full", title.expr = expression(GA[4]~"trap configurations (120 traps)"))

#GA5
GA5.120 <- bind_rows(ga51b, ga52b, ga53b, ga54b)

GA5.120 <- GA5.120 %>%
  mutate(
    design_label = case_when(
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[5]~S[1]",
      "GA[5]~S[2]",
      "GA[5]~Avg",
      "GA[5]~Both"
    ))
  )

GA5.120.plot <- plot.design(GA5.120, msk, view = "full", title.expr = expression(GA[5]~"trap configurations (120 traps)"))

ggsave("15FoldScen/figures/GA4120.pdf", plot = GA4.120.plot, width = 9, height = 6)
ggsave("15FoldScen/figures/GA5120.pdf", plot = GA5.120.plot, width = 9, height = 6)

#and both sets in one plot
#current fn no lomnger works for this so need to adjust if want GA4 and GA5 in 2 * 4 grid
opt.120 <- bind_rows(ga41b, ga42b, ga43b, ga44b,
                    ga51b, ga52b, ga53b, ga54b)

opt.120 <- opt.120 %>%
  mutate(
    design_label = case_when(
      design == "GA4 S1" ~ "GA[4]~S[1]",
      design == "GA4 S2" ~ "GA[4]~S[2]",
      design == "GA4 Avg" ~ "GA[4]~Av",
      design == "GA4 Both" ~ "GA[4]~Both",
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Av",
      design == "GA5 Both" ~ "GA[5]~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[4]~S[1]",
      "GA[5]~S[1]",
      "GA[4]~S[2]",
      "GA[5]~S[2]",
      "GA[4]~Av",
      "GA[5]~Av",
      "GA[4]~Both",
      "GA[5]~Both"
    ))
  )

opt.120.plot <- plot.design(opt.120, msk, view = "full", 
                            ndim1 = 2, ndim2 = 4, title.expr = expression(GA[4]~"and"~GA[5]~"trap configurations (120 traps)"))

ggsave("15FoldScen/figures/Opt120.pdf", plot = opt.120.plot, width = 7, height = 9.5)

##################################################
#put GA5 and 2 stage together
GA2.120 <- bind_rows(ga51b, ga52b, ga53b, ga54b,
                    ga2b)

GA2.120 <- GA2.120 %>%
  mutate(
    design_label = case_when(
      design == "GA5 S1" ~ "GA[5]~S[1]",
      design == "GA5 S2" ~ "GA[5]~S[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[5]~S[1]",
      "GA[5]~S[2]",
      "GA[5]~Avg",
      "GA[5]~Both",
      "Two~Stage"
    ))
  )

#now adding dummy panels to get proper alignment
sys.pad.120 <- bind_rows(
  sys.120,
  tibble(
    x = 0, y = 0,
    design = "dummy",
    design_label = "dummy"
  )
)

levels_all_sys <- c(levels(sys.120$design_label), "dummy")
sys.pad.120$design_label <- factor(sys.pad.120$design_label, levels = levels_all_sys)

ga4.pad.120 <- bind_rows(
  GA4.120,
  tibble(
    x = 0, y = 0,
    design = "dummy",
    design_label = "dummy"
  )
)

levels_all_ga4 <- c(levels(GA4.120$design_label), "dummy")
ga4.pad.120$design_label <- factor(ga4.pad.120$design_label, levels = levels_all_ga4)

#create trap locs with slightly larger extent to prevent LW spilling over
obj.lw <- create.extent(sigma = 3000, buff.factor = 2.900, res = 200)
trap.locs.lw <- obj.lw [[2]]

#use slightly larger trap loc extent to prevent LW traps spilling out of panel 
p.GA2.120 <- plot.design(GA2.120, trap.locs.lw, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL)
p.syspad.120 <- plot.design(sys.pad.120, trap.locs.lw, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL, levels_all = levels_all_sys)
p.GA4pad.120 <- plot.design(ga4.pad.120, trap.locs.lw, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL, levels_all = levels_all_ga4)

all.120.plot <- Combine.layoutplots(p.syspad.120, p.GA4pad.120,p.GA2.120)

setwd("~/Git/SCR-Design")
ggsave("15FoldScen/figures/All120.pdf", plot = all.120.plot, width = 270, height = 190, units = "mm")

