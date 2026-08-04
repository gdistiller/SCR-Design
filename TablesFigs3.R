#For tables, this script uses kableextra to construct tables in R
#then the latex code behind the table is saved to a .tex file
#lastly the code is included in Overleaf using: \input{tables/table1.tex}
#July 2026, updated with new OS and the new LW 120 
#uses AllResults4.RData
#updated Aug to use objkects for labels, colour scales etc

library(kableExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(terra)

#####################
# table of parameter estimates per strata
#####################
true.vals <- matrix(c(
  0.05,  0.05/15,    # D.x, D.y
  2,  2/15,  # L0.x, L0.y
  200, 3000 # Sig.x, Sig.y
), nrow = 3, byrow = TRUE)


rownames(true.vals) <- c("$D$", "$\\lambda_0$", "$\\sigma$")
colnames(true.vals) <- c("Group 1", "Group 2")

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

#####################
# % excluded
#####################

load("15FoldScen/AllResults4.RData")
AllResults40 <- Allresults_4[[1]]

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

df_excl <- AllResults40%>%
  group_by(sim, design, stratum) %>%
  summarise(
    bad_sim = any(is.na(estimate)),
    .groups = "drop"
  ) %>%
  group_by(design, stratum) %>%
  summarise(
    n_bad   = sum(bad_sim),
    n_sims  = n(),
    prop_bad = n_bad / n_sims,
    .groups = "drop"
  )

excl40.plot <- ggplot(df_excl,
       aes(x = design,
           y = prop_bad,
           fill = stratum)) +
  
  # ✅ bars (side-by-side within each design)
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7
  ) +
  
  # ✅ nice colours for strata
  scale_fill_manual(
    name = "Group",
    values = c(
      "S1" = "#1b9e77",
      "S2" = "#d95f02"
    ),
    labels = c(
      "S1" = "Group 1",
      "S2" = "Group 2"
    )
  ) +
  
  # ✅ show % on axis
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  
  # ✅ theme
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "top",
    strip.text = element_text(size = 12)
  ) +
  
  # ✅ labels
  labs(
    y = "Proportion excluded",
    x = "Design",
    title = "Simulation exclusion rates"
  )

ggsave("15FoldScen/figures/Excl40.pdf", plot = excl40.plot, width = 85, height = 85, units = "mm")

AllResults120 <- Allresults_4[[2]]

#confirm nothing dropped with 120 traps
#control order of design factor
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

df_excl_120 <- AllResults120%>%
  group_by(sim, design, stratum) %>%
  summarise(
    bad_sim = any(is.na(estimate)),
    .groups = "drop"
  ) %>%
  group_by(design, stratum) %>%
  summarise(
    n_bad   = sum(bad_sim),
    n_sims  = n(),
    prop_bad = n_bad / n_sims,
    .groups = "drop"
  )

df_excl_120$prop_bad

#################################################################
#figs created in ggplot and saved
#needs make.summary, Metric.plot and Combine.plots fns

#create summaries
summary.40 <- make.summary(AllResults40)
summary.120 <- make.summary(AllResults120)
summary.40$traps  <- "40"
summary.120$traps <- "120"
summary.combined <- dplyr::bind_rows(summary.40, summary.120)
summary.combined$traps <- factor(summary.combined$traps, levels = c("40", "120"))

#setup labels, scales etc
design_cols <- c(
  "Grid"     = "#0072B2",
  "Cluster"  = "#D55E00",
  "Lacework" = "#009E73",
  "min(n,r)"      = "#CC79A7",
  "En2"      = "#E69F00",
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

#with traps faceted
D.RB <- Metric.plot(summary.combined,
                     plot.title = NULL,
                     param_select = "Density",
                     metric = "RB",
                     facet_traps = TRUE,
                     ylims=c(-0.1,1))

#Rel SE plots (model-based)
D.RSE <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "Density",
                      metric = "RSE",
                      facet_traps = TRUE,
                      ylims=c(0,1))

#coverage plots
D.cov <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "Density",
                      metric = "coverage",
                      facet_traps = TRUE,
                      ylims=c(0.75,1))

#combine in two different ways
plot.D1 <- Combine.plots(D.RB, D.RSE, D.cov,
              global_title = "Performance for Density"
)

ggsave("15FoldScen/figures/D1.pdf", plot = plot.D1, width = 7, height = 9)

#2nd one old from when was trying diff plots
#ggsave("15FoldScen/figures/D2.pdf", plot = plot.D2, width = 7, height = 9)

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
                     ylims=c(0.75,1))

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
                     ylims=c(-0.1,1))

#Rel SE plots (model-based)
sig.RSE <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "sigma",
                      metric = "RSE",
                      facet_traps = TRUE,
                      ylims=c(0,1))

#coverage plots
sig.cov <- Metric.plot(summary.combined,
                      plot.title = NULL,
                      param_select = "sigma",
                      metric = "coverage",
                      facet_traps = TRUE,
                      ylims=c(0.75,1))

#combine in two different ways
plot.sig1 <- Combine.plots(sig.RB, sig.RSE, sig.cov,
                          global_title = expression("Performance of "~ sigma)
)
ggsave("15FoldScen/figures/Sig1.pdf", plot = plot.sig1, width = 7, height = 9)


###################################
#a plot of different realisations, for 40 traps
#objs with various proposed designs pasted into dir below
setwd("~/Git/SCR-Design/15FoldScen/Cluster/ProposedDesigns")
setwd("~/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns")

load("SCRObjs.RData")

load("GridDesigns40.RData")
load("LWdesigns.RData")
load("GADesigns40.RData")
load("GA2StageDesigns.RData") 

#extract coordinates and put all in one df

# extract designs, exclude the OS grid as only need to show one
sys1 <- as.data.frame(grid.designs.40$`800 m`[[2]]) %>%
  mutate(design = "Grid 800")

sys2 <- as.data.frame(grid.designs.40$`Cluster (opt)`[[1]]) %>%
  mutate(design = "Cluster (OS)")

sys3 <- as.data.frame(grid.designs.40$`Cluster (2 sig)`[[1]]) %>%
  mutate(design = "Cluster (2σ)")

sys4 <- as.data.frame(lwlist$`40 traps`) %>%
  mutate(design = "Lacework")

ga41 <- as.data.frame(GA.designs.40$G4S1[[1]]) %>%
  mutate(design = "min(n,r)-G1")

ga42 <- as.data.frame(GA.designs.40$G4S2[[1]]) %>%
  mutate(design = "min(n,r)-G2")

ga43 <- as.data.frame(GA.designs.40$G4Avg[[1]]) %>%
  mutate(design = "min(n,r)-A")

ga44 <- as.data.frame(GA.designs.40$G4Both[[1]]) %>%
  mutate(design = "min(n,r)-B")

ga51 <- as.data.frame(GA.designs.40$G5S1[[1]]) %>%
  mutate(design = "En2-G1")

ga52 <- as.data.frame(GA.designs.40$G5S2[[1]]) %>%
  mutate(design = "En2-G2")

ga53 <- as.data.frame(GA.designs.40$G5Avg[[1]]) %>%
  mutate(design = "En2-A")

ga54 <- as.data.frame(GA.designs.40$G5Both[[1]]) %>%
  mutate(design = "En2-B")

ga2 <- as.data.frame(GA2StageDesigns$`40 traps`[[1]]) %>%
  mutate(design = "Two Stage")

#first deal with systematic designs
sys.40 <- bind_rows(sys1, sys2, sys3, sys4)

#set labels to parse properly
sys.40 <- sys.40 %>%
  mutate(
    design_label = case_when(
      design == "Grid (OS)" ~ "Grid~'(OS)'",
      design == "Grid 800" ~ "Grid~800",
      design == "Cluster (OS)" ~ "Cluster~'(OS)'",
      design == "Cluster (2σ)" ~ "Cluster~(2*sigma)",
      design == "Lacework" ~ "Lacework"
    ),
    design_label = factor(design_label, levels = c(
      "Grid~'(OS)'",
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
setwd("~/Documents/Git/SCRDesign_fresh")
ggsave("15FoldScen/figures/sys40.pdf", plot = sys40.plot, width = 9, height = 6)

################################################################################
#now optimised designs plotted on full extent
#first just GA4 or GA5
GA4.40 <- bind_rows(ga41, ga42, ga43, ga44)

GA4.40 <- GA4.40 %>%
  mutate(
    design_label = case_when(
      design == "min(n,r)-G1" ~ "'min(n,r)'~G[1]",
      design == "min(n,r)-G2" ~ "'min(n,r)'~G[2]",
      design == "min(n,r)-A" ~ "'min(n,r)'~Avg",
      design == "min(n,r)-B" ~ "'min(n,r)'~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "'min(n,r)'~G[1]",
      "'min(n,r)'~G[2]",
      "'min(n,r)'~Avg",
      "'min(n,r)'~Both")
      )
  )

GA4.40.plot <- plot.design(GA4.40, msk,  view = "full", ndim2 = 4, title.expr = expression(GA[4]~"trap configurations (40 traps)"))

#GA5
GA5.40 <- bind_rows(ga51, ga52, ga53, ga54)

GA5.40 <- GA5.40 %>%
  mutate(
    design_label = case_when(
      design == "En2-G1" ~ "En2~G[1]",
      design == "En2-G2" ~ "En2~G[2]",
      design == "En2-A" ~ "En2~Avg",
      design == "En2-B" ~ "En2~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "En2~G[1]",
      "En2~G[2]",
      "En2~Avg",
      "En2~Both"
    ))
  )

GA5.40.plot <- plot.design(GA5.40, msk, view = "full", ndim2 = 4, title.expr = expression(GA[5]~"trap configurations (40 traps)"))

ggsave("15FoldScen/figures/GA440.pdf", plot = GA4.40.plot, width = 9, height = 6)
ggsave("15FoldScen/figures/GA540.pdf", plot = GA5.40.plot, width = 9, height = 6)

##################################################################################
#put GA5 and 2 stage together, redundant now
GA2.40 <- bind_rows(ga51, ga52, ga53, ga54,
                     ga2)

GA2.40 <- GA2.40 %>%
  mutate(
    design_label = case_when(
      design == "GA5 G1" ~ "GA[5]~G[1]",
      design == "GA5 G2" ~ "GA[5]~G[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[5]~G[1]",
      "GA[5]~G[2]",
      "GA[5]~Avg",
      "GA[5]~Both",
      "Two~Stage"
    ))
  )

#put GA4 and 2 stage together
GA2b.40 <- bind_rows(ga41, ga42, ga43, ga44,
                    ga2)

GA2b.40 <- GA2b.40 %>%
  mutate(
    design_label = case_when(
      design == "min(n,r)-G1" ~ "'min(n,r)'~G[1]",
      design == "min(n,r)-G2" ~ "'min(n,r)'~G[2]",
      design == "min(n,r)-A" ~ "'min(n,r)'~Avg",
      design == "min(n,r)-B" ~ "'min(n,r)'~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "'min(n,r)'~G[1]",
      "'min(n,r)'~G[2]",
      "'min(n,r)'~Avg",
      "'min(n,r)'~Both",
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

ga5.pad.40 <- bind_rows(
  GA5.40,
  tibble(
    x = 0, y = 0,
    design = "dummy",
    design_label = "dummy"
  )
)

levels_all_ga5 <- c(levels(GA5.40$design_label), "dummy")
ga5.pad.40$design_label <- factor(ga5.pad.40$design_label, levels = levels_all_ga5)

#p.GA2.40 <- plot.design(GA2.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, 
#                        point.size = 0.5, title.expr = NULL)
p.GA2b.40 <- plot.design(GA2b.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, 
                        point.size = 0.5, title.expr = NULL)
p.syspad.40 <- plot.design(sys.pad.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, title.expr = NULL, 
                           point.size = 0.5, buffer.prop = 0.05, levels_all = levels_all_sys)
#p.GA4pad.40 <- plot.design(ga4.pad.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, 
#                           point.size = 0.5, title.expr = NULL, levels_all = levels_all_ga4)
p.GA5pad.40 <- plot.design(ga5.pad.40, msk.red, view = "full", ndim1 = 1, ndim2 = 5, 
                           point.size = 0.5, title.expr = NULL, levels_all = levels_all_ga5)

#all.40.plot <- Combine.layoutplots(p.syspad.40, p.GA4pad.40,p.GA2.40)
all.40.plotb <- Combine.layoutplots(p.syspad.40, p.GA2b.40, p.GA5pad.40)

setwd("~/Git/SCR-Design")
#ggsave("15FoldScen/figures/All40.pdf", plot = all.40.plot, width = 270, height = 190, units = "mm")
ggsave("15FoldScen/figures/All40.pdf", plot = all.40.plotb, width = 270, height = 190, units = "mm")

###############
#for 120 traps
#version xxb is from ngen = 500
###############
setwd("~/Documents/Git/SCRDesign_fresh/15FoldScen/Cluster/Sims/120Traps")
setwd("~/Git/SCR-Design/15FoldScen/Cluster/Sims/120Traps")

load("GridDesigns120.RData")
load("GADesigns120.RData")

setwd("~/Documents/Git/SCRDesign_fresh/15FoldScen/Cluster/Sims")
setwd("~/Git/SCR-Design/15FoldScen/Cluster/Sims")
load("LWdesigns.RData")
load("GA2StageDesignsb.RData") 

#extract coordinates and put all in one df

# extract designs, again exclude the OS grid
sys1b <- as.data.frame(grid.designs.120$`800 m`[[2]]) %>%
  mutate(design = "Grid 800")

sys2b <- as.data.frame(grid.designs.120$`Cluster (opt)`[[1]]) %>%
  mutate(design = "Cluster (OS)")

sys3b <- as.data.frame(grid.designs.120$`Cluster (2 sig)`[[1]]) %>%
  mutate(design = "Cluster (2σ)")

sys4b <- as.data.frame(lwlist$`120 traps`) %>%
  mutate(design = "Lacework")

sys5b <- as.data.frame(lwlist$`120 trapsB`) %>%
  mutate(design = "Lacework (F)")

ga41b <- as.data.frame(GA.designs.120$G4S1[[1]]) %>%
  mutate(design = "min(n,r)-G1")

ga42b <- as.data.frame(GA.designs.120$G4S2[[1]]) %>%
  mutate(design = "min(n,r)-G2")

ga43b <- as.data.frame(GA.designs.120$G4Avg[[1]]) %>%
  mutate(design = "min(n,r)-A")

ga44b <- as.data.frame(GA.designs.120$G4Both[[1]]) %>%
  mutate(design = "min(n,r)-B")

ga51b <- as.data.frame(GA.designs.120$G5S1[[1]]) %>%
  mutate(design = "En2-G1")

ga52b <- as.data.frame(GA.designs.120$G5S2[[1]]) %>%
  mutate(design = "En2-G2")

ga53b <- as.data.frame(GA.designs.120$G5Avg[[1]]) %>%
  mutate(design = "En2-A")

ga54b <- as.data.frame(GA.designs.120$G5Both[[1]]) %>%
  mutate(design = "En2-B")

ga2b <- as.data.frame(GA2StageDesigns$`120 traps`[[1]]) %>%
  mutate(design = "Two Stage")

#first deal with systematic designs
sys.120 <- bind_rows(sys1b, sys2b, sys3b, sys4b, sys5b)

#set labels to parse properly
sys.120 <- sys.120 %>%
  mutate(
    design_label = case_when(
      design == "Grid 800" ~ "Grid~800",
      design == "Cluster (OS)" ~ "Cluster~'(OS)'",
      design == "Cluster (2σ)" ~ "Cluster~(2*sigma)",
      design == "Lacework" ~ "Lacework",
      design == "Lacework (F)" ~ "Lacework (F)"
    ),
    design_label = factor(design_label, levels = c(
      "Grid~800",
      "Cluster~'(OS)'",
      "Cluster~(2*sigma)",
      "Lacework",
      "Lacework (F)"
    ))
  )

msk <- SCR.objs$Full[[1]]

#use fn
#move lacework coords up a few hundred meters
sys.120$y[sys.120$design=="Lacework"] <- sys.120$y[sys.120$design=="Lacework"] + 1000
sys120.plot <- plot.design(sys.120, msk, view = "full", title = "Systematic trap configurations (120 traps)")

setwd("~/Documents/Git/SCRDesign_fresh")
setwd("~/Git/SCR-Design")
ggsave("15FoldScen/figures/sys120.pdf", plot = sys120.plot, width = 9, height = 6)

#######################################################################################
#now optimised designs plotted on full extent
#first just GA4 or GA5
GA4.120 <- bind_rows(ga41b, ga42b, ga43b, ga44b)

GA4.120 <- GA4.120 %>%
  mutate(
    design_label = case_when(
      design == "min(n,r)-G1" ~ "'min(n,r)'~G[1]",
      design == "min(n,r)-G2" ~ "'min(n,r)'~G[2]",
      design == "min(n,r)-A" ~ "'min(n,r)'~Avg",
      design == "min(n,r)-B" ~ "'min(n,r)'~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "'min(n,r)'~G[1]",
      "'min(n,r)'~G[2]",
      "'min(n,r)'~Avg",
      "'min(n,r)'~Both")
    )
  )

GA4.120.plot <- plot.design(GA4.120, msk, view = "full", title.expr = expression(GA[4]~"trap configurations (120 traps)"))

#GA5
GA5.120 <- bind_rows(ga51b, ga52b, ga53b, ga54b)

GA5.120 <- GA5.120 %>%
  mutate(
    design_label = case_when(
      design == "En2-G1" ~ "En2~G[1]",
      design == "En2-G2" ~ "En2~G[2]",
      design == "En2-A" ~ "En2~Avg",
      design == "En2-B" ~ "En2~Both",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "En2~G[1]",
      "En2~G[2]",
      "En2~Avg",
      "En2~Both"
    ))
  )

GA5.120.plot <- plot.design(GA5.120, msk, view = "full", title.expr = expression(GA[5]~"trap configurations (120 traps)"))

ggsave("15FoldScen/figures/GA4120.pdf", plot = GA4.120.plot, width = 9, height = 6)
ggsave("15FoldScen/figures/GA5120.pdf", plot = GA5.120.plot, width = 9, height = 6)

##################################################################################
#put GA5 and 2 stage together, redundant
GA2.120 <- bind_rows(ga51b, ga52b, ga53b, ga54b,
                    ga2b)

GA2.120 <- GA2.120 %>%
  mutate(
    design_label = case_when(
      design == "GA5 G1" ~ "GA[5]~G[1]",
      design == "GA5 G2" ~ "GA[5]~G[2]",
      design == "GA5 Avg" ~ "GA[5]~Avg",
      design == "GA5 Both" ~ "GA[5]~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "GA[5]~G[1]",
      "GA[5]~G[2]",
      "GA[5]~Avg",
      "GA[5]~Both",
      "Two~Stage"
    ))
  )

#put GA4 and 2 stage together, optimisation of stage 2 is GA4
GA2b.120 <- bind_rows(ga41b, ga42b, ga43b, ga44b,
                     ga2b)

GA2b.120 <- GA2b.120 %>%
  mutate(
    design_label = case_when(
      design == "min(n,r)-G1" ~ "'min(n,r)'~G[1]",
      design == "min(n,r)-G2" ~ "'min(n,r)'~G[2]",
      design == "min(n,r)-A" ~ "'min(n,r)'~Avg",
      design == "min(n,r)-B" ~ "'min(n,r)'~Both",
      design == "Two Stage" ~ "Two~Stage",
      TRUE ~ design
    ),
    design_label = factor(design_label, levels = c(
      "'min(n,r)'~G[1]",
      "'min(n,r)'~G[2]",
      "'min(n,r)'~Avg",
      "'min(n,r)'~Both",
      "Two~Stage"
    ))
  )

#now adding dummy panels to get proper alignment
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

ga5.pad.120 <- bind_rows(
  GA5.120,
  tibble(
    x = 0, y = 0,
    design = "dummy",
    design_label = "dummy"
  )
)

levels_all_ga5 <- c(levels(GA5.120$design_label), "dummy")
ga5.pad.120$design_label <- factor(ga5.pad.120$design_label, levels = levels_all_ga5)

#decided to shift coords for LW traps up rather than adjust buffer
#create trap locs to plot without buffer
obj <- create.extent(sigma = 3000, buff.factor = 3, res = 200)
trap.locs <- obj [[2]]

###########################################################
#have put 2 stage with GA4, still displaying in middle row for now

#p.GA2.120 <- plot.design(GA2.120, trap.locs, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL)
p.GA2b.120 <- plot.design(GA2b.120, trap.locs, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL)
#p.GA4pad.120 <- plot.design(ga4.pad.120, trap.locs, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL, levels_all = levels_all_ga4)
p.GA5pad.120 <- plot.design(ga5.pad.120, trap.locs, view = "full", ndim1 = 1, ndim2 = 5, point.size = 0.5, title.expr = NULL, levels_all = levels_all_ga5)

#all.120.plot <- Combine.layoutplots(sys120.plot, p.GA4pad.120,p.GA2.120)
all.120.plotb <- Combine.layoutplots(sys120.plot, p.GA2b.120, p.GA5pad.120)

setwd("~/Git/SCR-Design")
#ggsave("15FoldScen/figures/All120.pdf", plot = all.120.plot, width = 270, height = 190, units = "mm")
ggsave("15FoldScen/figures/All120.pdf", plot = all.120.plotb, width = 270, height = 190, units = "mm")

