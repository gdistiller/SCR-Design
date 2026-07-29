#July 2026, process sim results for 15 fold scen, 500 reps
#filters out reps where se/est > 100, uses mag factor of ten fold
#excludes GA both (max-min) and other spacings not used
#new OS for cluster and grid 120
#adds new LW 120 design

library(secr)
library(secrdesign)
library(dplyr)
library(tidyr)

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

####################################################
#set some global values for plotting and filtering
mag.factor <- 10
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

#600 m spacing (OS)
load("15FoldScen/Cluster/Sims/40Traps/GridseResults40.RData")

Grid40.600.red1 <- find.rogue(Grid.600.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid40.600.red2 <- find.rogue(Grid.600.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid40.600.clean <- merge(Grid40.600.red1, Grid40.600.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid600.40 <- construct_df(Grid40.600.clean, design_name = "Grid (OS)", true_vals = truth_df)

#800 m spacing

load("15FoldScen/Cluster/Sims/40Traps/GridsdResults40.RData")

Grid40.800.red1 <- find.rogue(Grid.800.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid40.800.red2 <- find.rogue(Grid.800.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid40.800.clean <- merge(Grid40.800.red1, Grid40.800.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid800.40 <- construct_df(Grid40.800.clean, design_name = "Grid 800m", true_vals = truth_df)

#cluster (OS spacing)
load("15FoldScen/Cluster/Sims/40Traps/OSclusters40.RData")

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

#########
#Lacework
#########

load("15FoldScen/Cluster/Sims/40Traps/LWResults40.RData")

LW40.red1 <- find.rogue(LW.40.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LW40.red2 <- find.rogue(LW.40.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
LW40.clean <- merge(LW40.red1, LW40.red2, by = "row", all.x = TRUE, all.y = TRUE)

LW.40 <- construct_df(LW40.clean, design_name = "Lacework", true_vals = truth_df)

##########
#Two Stage
##########

load("15FoldScen/Cluster/Sims/40Traps/TwoStageResults40.RData")

TwoStage40.red1 <- find.rogue(TwoStage.40.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
TwoStage40.red2 <- find.rogue(TwoStage.40.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
TwoStage40.clean <- merge(TwoStage40.red1, TwoStage40.red2, by = "row", all.x = TRUE, all.y = TRUE)

TwoStage.40 <- construct_df(TwoStage40.clean, design_name = "Two Stage", true_vals = truth_df)

all_results_40_4 <- bind_rows(
  Grid600.40,
  Grid800.40,
  ClusterOS.40,
  Cluster2Sig.40,
  LW.40,
  GA4.S1.40,
  GA4.S2.40,
  GA4.Avg.40,
  GA4.Both.40,
  GA5.S1.40,
  GA5.S2.40,
  GA5.Avg.40,
  GA5.Both.40,
  TwoStage.40
  )

#####################################

#################
#120 traps
#################

###########################################################
#Grid designs (incl grids with OS 2G, and clustered grids)
###########################################################

#1300 m spacing (OS)
#ran as 5 jobs, need to collate

load("15FoldScen/Cluster/Sims/120Traps/OSGridsResults120a.RData")
load("15FoldScen/Cluster/Sims/120Traps/OSGridsResults120b.RData")
load("15FoldScen/Cluster/Sims/120Traps/OSGridsResults120c.RData")
load("15FoldScen/Cluster/Sims/120Traps/OSGridsResults120d.RData")
load("15FoldScen/Cluster/Sims/120Traps/OSGridsResults120e.RData")


Grid.1300.results <- rbind(
  Grid.1300.results.a,
  Grid.1300.results.b,
  Grid.1300.results.c,
  Grid.1300.results.d,
  Grid.1300.results.e
)

Grid120.1300.red1 <- find.rogue(Grid.1300.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid120.1300.red2 <- find.rogue(Grid.1300.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid120.1300.clean <- merge(Grid120.1300.red1, Grid120.1300.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid1300.120 <- construct_df(Grid120.1300.clean, design_name = "Grid (OS)", true_vals = truth_df)

#800 m spacing

load("15FoldScen/Cluster/Sims/120Traps/Grids800Results120.RData")

Grid120.800.red1 <- find.rogue(Grid.800.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid120.800.red2 <- find.rogue(Grid.800.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
Grid120.800.clean <- merge(Grid120.800.red1, Grid120.800.red2, by = "row", all.x = TRUE, all.y = TRUE)

Grid800.120 <- construct_df(Grid120.800.clean, design_name = "Grid 800m", true_vals = truth_df)

#cluster (os spacing)
load("15FoldScen/Cluster/Sims/120Traps/OSclusters120.RData")

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

#########
#Lacework
#########

load("15FoldScen/Cluster/Sims/120Traps/LWResults120.RData")

LW120.red1 <- find.rogue(LW.120.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LW120.red2 <- find.rogue(LW.120.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
LW120.clean <- merge(LW120.red1, LW120.red2, by = "row", all.x = TRUE, all.y = TRUE)

LW.120 <- construct_df(LW120.clean, design_name = "Lacework", true_vals = truth_df)

load("15FoldScen/Cluster/Sims/120Traps/LWBResults120.RData")

LWB120.red1 <- find.rogue(LWB.120.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LWB120.red2 <- find.rogue(LWB.120.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
LWB120.clean <- merge(LWB120.red1, LWB120.red2, by = "row", all.x = TRUE, all.y = TRUE)

LWB.120 <- construct_df(LWB120.clean, design_name = "Lacework (F)", true_vals = truth_df)

##########
#Two Stage
##########

load("15FoldScen/Cluster/Sims/120Traps/TwoStageResults120b.RData")

TwoStage120.red1 <- find.rogue(TwoStage.120.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
TwoStage120.red2 <- find.rogue(TwoStage.120.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))
TwoStage120.clean <- merge(TwoStage120.red1, TwoStage120.red2, by = "row", all.x = TRUE, all.y = TRUE)

TwoStage.120 <- construct_df(TwoStage120.clean, design_name = "Two Stage", true_vals = truth_df)


################################################################

all_results_120_4 <- bind_rows(
  Grid1300.120,
  Grid800.120,
  ClusterOS.120,
  Cluster2Sig.120,
  LW.120,
  LWB.120,
  GA4.S1.120,
  GA4.S2.120,
  GA4.Avg.120,
  GA4.Both.120,
  GA5.S1.120,
  GA5.S2.120,
  GA5.Avg.120,
  GA5.Both.120,
  TwoStage.120)

Allresults_4 <- list("40 traps" = all_results_40_4, "120 traps" = all_results_120_4)
save(Allresults_4, file = "15FoldScen/AllResults4.RData")

########################################################################
