#June 2026, process sim results for 15 fold scen, 500 reps
#filters out reps where se/est > 100
#and this version increases the mag factor to ten fold
#also excludes GA both (max-min)

library(secr)
library(secrdesign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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

TwoStage.40 <- construct_df(TwoStage40.clean, design_name = "TwoStage", true_vals = truth_df)


all_results_40_3 <- bind_rows(
  Grid800.40,
  Grid700.40,
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
  TwoStage
  )

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

################################################################

all_results_120_3 <- bind_rows(
  Grid800.120,
  Grid700.120,
  ClusterOS.120,
  Cluster2Sig.120,
  LW.120,
  GA4.S1.120,
  GA4.S2.120,
  GA4.Avg.120,
  GA4.Both.120,
  GA5.S1.120,
  GA5.S2.120,
  GA5.Avg.120,
  GA5.Both.120)

Allresults_3 <- list("40 traps" = all_results_40_3, "120 traps" = all_results_120_3)
save(Allresults_3, file = "15FoldScen/AllResults3.RData")

########################################################################
#trying different plots
load("15FoldScen/AllResults3.RData")
AllResults40 <- Allresults_3[[1]]
AllResults120 <- Allresults_3[[2]]
AllResults40$traps  <- "40"
AllResults120$traps <- "120"
AllResults.combined <- bind_rows(AllResults40, AllResults120)
AllResults.combined$traps <- factor(AllResults.combined$traps, levels = c("40", "120"))

#code to generate bar plot of exclusions, first constructs df
df_excl <- AllResults.combined %>%
  filter(traps == "40") %>%   # <- keep only 40 trap
  group_by(sim, design, stratum, traps) %>%
  summarise(
    bad_sim = any(is.na(estimate)),
    .groups = "drop"
  ) %>%
  group_by(design, stratum, traps) %>%
  summarise(
    n_bad   = sum(bad_sim),
    n_sims  = n(),
    prop_bad = n_bad / n_sims,
    .groups = "drop"
  )

ggplot(df_excl,
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
    name = "Stratum",
    values = c(
      "S1" = "#1b9e77",
      "S2" = "#d95f02"
    ),
    labels = c(
      "S1" = "Stratum 1",
      "S2" = "Stratum 2"
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

#now plots of RB, Rel SE (model based) and coverage
#Done per parameter and patched together

#create summaries
summary.40 <- make.summary(AllResults40)
summary.120 <- make.summary(AllResults120)
summary.40$traps  <- "40"
summary.120$traps <- "120"
summary.combined <- dplyr::bind_rows(summary.40, summary.120)
summary.combined$traps <- factor(summary.combined$traps, levels = c("40", "120"))

plot_metric(summary.combined,
            yvar = "RB",
            sevar = "RB_se",
            ylab = "Relative Bias",
            yline = 0,
            plot.title = "Relative Bias",
            ylims = c(-0.1, 0.6))

#RB plots per parameter
D.RB <- Metric.plot(summary.combined,
                  plot.title = "Relative Bias (Density)",
                  param_select = "Density",
                  metric = "RB",
                  ylims = c(-0.2, 1))

L0.RB <- Metric.plot(summary.combined,
                    plot.title = expression("Relative Bias (" * lambda[0] * ")"),
                    param_select = "lambda[0]",
                    metric = "RB",
                    ylims = c(-0.2, 0.2))  

Sig.RB <- Metric.plot(summary.combined,
                     plot.title = expression("Relative Bias (" * sigma * ")"),
                     param_select = "sigma",
                     metric = "RB",
                     ylims = c(-0.2, 1))  

#Rel SE plots (model-based)
D.RSE <- Metric.plot(summary.combined,
                    plot.title = "Relative SE (Density)",
                    param_select = "Density",
                    metric = "RSE",
                    ylims = c(-0.2, 1))

L0.RSE <- Metric.plot(summary.combined,
                     plot.title = expression("Relative SE (" * lambda[0] * ")"),
                     param_select = "lambda[0]",
                     metric = "RSE",
                     ylims = c(-0.2, 0.2))  

Sig.RSE <- Metric.plot(summary.combined,
                      plot.title = expression("Relative SE (" * sigma * ")"),
                      param_select = "sigma",
                      metric = "RSE",
                      ylims = c(-0.2, 1))  

#coverage plots
D.cov <- Metric.plot(summary.combined,
                     plot.title = "Coverage (Density)",
                     param_select = "Density",
                     metric = "coverage",
                     ylims = c(0.7, 1))

L0.cov <- Metric.plot(summary.combined,
                      plot.title = expression("Coverage (" * lambda[0] * ")"),
                      param_select = "lambda[0]",
                      metric = "coverage",
                      ylims = c(0.7, 1))  

Sig.cov <- Metric.plot(summary.combined,
                       plot.title = expression("Coverage (" * sigma * ")"),
                       param_select = "sigma",
                       metric = "coverage",
                       ylims = c(0.7, 1))  

#try new version
Metric.plot(summary.combined,
            plot.title = "RB (Density)",
            param_select = "Density",
            metric = "RB",
            facet_traps = TRUE)


D.RB <- Metric.plot(summary.combined,
            plot.title = "RB (Density)",
            param_select = "Density",
            metric = "RB",
            facet_traps = FALSE,
            ylims=c(-0.1,1))

D.RSE <- Metric.plot(summary.combined,
                    plot.title = "RSE (Density)",
                    param_select = "Density",
                    metric = "RSE",
                    facet_traps = FALSE,
                    ylims=c(0,1))

D.cov <- Metric.plot(summary.combined,
                     plot.title = "Coverage (Density)",
                     param_select = "Density",
                     metric = "coverage",
                     facet_traps = FALSE,
                     ylims=c(0.75,1))




p_combined <-
  (D.RB / D.RSE / D.cov) +
  plot_layout(guides = "collect", heights = c(1, 1, 1)) &
  theme(
    legend.position = "top",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    
    axis.text.x = element_text(size = 7, angle = 70, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 9),
    
    strip.text = element_text(size = 9),
    
    panel.spacing = unit(0.3, "lines"),
    plot.margin = margin(2, 2, 2, 2)
  )

# ✅ Add global title + panel labels
p_combined <- p_combined +
  plot_annotation(
    title = "test",
    tag_levels = "A"
  )

p_combined 
