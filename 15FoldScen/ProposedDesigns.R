#28 March 2025, run again in November to get 500 reps
#getting proposed designs for new 15 fold diff scen
#this code generates the grid designs, the GA designs are done on the cluster
#cluster files are all combined into a GADesigns file at the end

library(secrdesign)
library(kofnGA)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

rm(list=ls())

#Set values for both strata
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
D1 = 0.05 ; D2 = D1 / DiffFactor  ; D <- c(D1, D2)
L01 <- 2 ; L02 <- L01/DiffFactor
nT <- 40 ; nreps = 500

################################
#45 km area
################################

#setup mask and trap locs
#create mask and potential traplocs
#using 15 * larger sigma for extent and 3 sigma for buffer
res.objs <- create.extent(sigma = 3000)
mask <- res.objs[[1]]
trap.locs <- res.objs[[2]]
traplocs.sf <- res.objs[[3]]

plot(mask, axes = T)
plot(trap.locs, add = T)

#create grid designs using 800 m
grid.800 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                             lambda0 = NULL, sigma.buff = NULL, grid.spacing = 800, 
                             criterion = 4, n.reps = nreps, grid = TRUE, nT = nT)

#1600
grid.1600<- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = mean(sigma), 
                              criterion = 4, n.reps = nreps, grid = TRUE, nT = nT)

grid.designs <- list("800 m" = grid.800, "Avg sigma" = grid.1600)

setwd("~/Library/CloudStorage/OneDrive-UniversityofCapeTown/Documents/Git/SCRDesign/15FoldScen")
save(grid.designs, file = "GridDesigns.RData")
save(res.objs, file = "SCRObjs.RData")

##########################################################
#collate GA designs
#now uses version _3 to get 500 reps
#initial run had timeout so was rerun to save separate objs for GA4 and GA5

load("15FoldScen/Cluster/GA4Strata1Pars3.RData")
load("15FoldScen/Cluster/GA5Strata1Pars3.RData")

load("15FoldScen/Cluster/GA4Strata2Pars3.RData")
load("15FoldScen/Cluster/GA5Strata2Pars3.RData")

load("15FoldScen/Cluster/GA4AvgStrataPars3.RData")
load("15FoldScen/Cluster/GA5AvgStrataPars3.RData")

#using both strata resulted in timeouts, hence a and b files for GA4 and GA5
#combine the two sets of 250 into a consolidated list of 500
#Both strata
load("15FoldScen/Cluster/GA4BothStrataPars3a.RData")
GA4.BothStrataPars.250reps <- GA4.BothStrataPars3.40.list
  
load("15FoldScen/Cluster/GA4BothStrataPars3b.RData")

GA4.BothStrataPars3.40.list <- c(GA4.BothStrataPars.250reps, GA4.BothStrataPars3.40.list)
rm(GA4.BothStrataPars.250reps)

load("15FoldScen/Cluster/GA5BothStrataPars3a.RData")
GA5.BothStrataPars.250reps <- GA5.BothStrataPars3.40.list

load("15FoldScen/Cluster/GA5BothStrataPars3b.RData")
GA5.BothStrataPars3.40.list <- c(GA5.BothStrataPars.250reps, GA5.BothStrataPars3.40.list)
rm(GA5.BothStrataPars.250reps)

#Both strata max-min
load("15FoldScen/Cluster/GA4BothStrataParsMax3a.RData")
GA4.BothStrataParsMax.250reps <- GA4.BothStrataParsMax3.40.list

load("15FoldScen/Cluster/GA4BothStrataParsMax3b.RData")
GA4.BothStrataParsMax3.40.list <- c(GA4.BothStrataParsMax.250reps, GA4.BothStrataParsMax3.40.list)
rm(GA4.BothStrataParsMax.250reps)

load("15FoldScen/Cluster/GA5BothStrataParsMax3a.RData")
GA5.BothStrataParsMax.250reps <- GA5.BothStrataParsMax3.40.list

load("15FoldScen/Cluster/GA5BothStrataParsMax3b.RData")
GA5.BothStrataParsMax3.40.list <- c(GA5.BothStrataParsMax.250reps, GA5.BothStrataParsMax3.40.list)
rm(GA5.BothStrataParsMax.250reps)

GA.designs <- list("G4S1" = GA4.Strata1Pars3.40.list, 
                   "G4S2" = GA4.Strata2Pars3.40.list,
                   "G5S1" = GA5.Strata1Pars3.40.list, 
                   "G5S2" = GA5.Strata2Pars3.40.list,
                   "G4Avg" = GA4.AvgStrataPars3.40.list, 
                   "G5Avg" = GA5.AvgStrataPars3.40.list,
                   "G4Both" = GA4.BothStrataPars3.40.list, 
                   "G5Both" = GA5.BothStrataPars3.40.list,
                   "G4BothMax" = GA4.BothStrataParsMax3.40.list, 
                   "G5BothMax" = GA5.BothStrataParsMax3.40.list)

save(GA.designs, file = "15FoldScen/GADesigns.RData")
