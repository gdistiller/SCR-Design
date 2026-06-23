#For tables, this script uses kableextra to construct tables in R
#then the latex code behind the table is saved to a .tex file
#lastly the code is included in oOverleaf using: \input{tables/table1.tex}

library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# table of parameter esrtimates per strata

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
  "Lacework" = c("Lacework",
                 "Two Stage" = c("2 Stage"))
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
                     plot.title = "RB (Density)",
                     param_select = "Density",
                     metric = "RB",
                     facet_traps = TRUE,
                     ylims=c(-0.1,1))

#Rel SE plots (model-based)
D.RSE2 <- Metric.plot(summary.combined,
                      plot.title = "RSE (Density)",
                      param_select = "Density",
                      metric = "RSE",
                      facet_traps = TRUE,
                      ylims=c(0,1))

#coverage plots
D.cov2 <- Metric.plot(summary.combined,
                      plot.title = "Coverage (Density)",
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



###################################
#a plot of different realisations, for 40 traps
#objs with various proposed designs pasted into dir below
setwd("~/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns")

load("SCRObjs.RData")

load("GridDesigns40.RData")
load("LWdesigns.RData")
load("GADesigns40.RData")
load("GA2StageDesigns.RData") 

#extract coordinates and put all in one df
# mask
mask_df <- as.data.frame(mask)  # should give x, y

# extract designs
sys1 <- as.data.frame(grid.designs.40$`800 m`[[1]]) %>%
  mutate(design = "Grid 800")

sys2 <- as.data.frame(grid.designs.40$`Cluster (opt)`[[1]]) %>%
  mutate(design = "Cluster (opt)")

sys3 <- as.data.frame(grid.designs.40$`Cluster (2 sig)`[[1]]) %>%
  mutate(design = "Cluster (2 sig)")

sys4 <- as.data.frame(lwlist$`40 traps`) %>%
  mutate(design = "Lacework")

ga1 <- as.data.frame(GA.designs.40$G4S1[[1]]) %>%
  mutate(design = "GA4 S1")

ga2 <- as.data.frame(GA.designs.40$G4S2[[1]]) %>%
  mutate(design = "GA4 S2")

ga3 <- as.data.frame(GA.designs.40$G4Avg[[1]]) %>%
  mutate(design = "GA4 Avg")

ga4 <- as.data.frame(GA.designs.40$G4Both[[1]]) %>%
  mutate(design = "GA4 Both")

ga5 <- as.data.frame(GA.designs.40$G5S1[[1]]) %>%
  mutate(design = "GA5 S1")

ga6 <- as.data.frame(GA.designs.40$G5S2[[1]]) %>%
  mutate(design = "GA5 S2")

ga7 <- as.data.frame(GA.designs.40$G5Avg[[1]]) %>%
  mutate(design = "GA5 Avg")

ga8 <- as.data.frame(GA.designs.40$G5Both[[1]]) %>%
  mutate(design = "GA5 Both")

ga9 <- as.data.frame(GA2StageDesigns$`40 traps`[[1]]) %>%
  mutate(design = "Two Stage")

design.40.df <- bind_rows(sys1, sys2, sys3, sys4,
                       ga1, ga2, ga3, ga4, ga5, ga6, ga7, ga8, ga9)


mask <- SCR.objs$Full[[1]]




