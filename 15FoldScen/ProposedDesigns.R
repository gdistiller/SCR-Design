#March 2026, incls OS and clustered designs
#this code generates the grid designs, the GA designs are done on the cluster
#cluster files are all combined into a GADesigns file at the end

library(secrdesign)
library(kofnGA)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

rm(list=ls())

source("functions2.R")

#Set values for both strata
sigma1 = 1500 ; sigma2 = 3000  ; sigma <- c(sigma1, sigma2)
D1 = 0.05 ; D2 = 0.01  ; D <- c(D1, D2)
L01 <- 2 ; L02 <- 0.5
nT <- 40 ; nreps = 500

################################
#GKE area
################################

#load masks created from shape files
#create mask and potential traplocs
#using 3 * larger sigma for buffer
load("GKE/GKEMasks.RData")

plot(GKEmsks$`GKE Mask`, axes = T)
plot(GKEmsks$`GKE Trap locs`, add = T)

#######################
#Otimal spacing grid
#create grid designs using optimal spacing, sum_min,  all_min critertia for two classes
#result does not seem sensitive to spacing of base grid
#######################

grid <- make.grid(7,7,1000, detector = "count")
spacing.sum.min <- optimalSpacing2G(D1 = D1, D2 = D2,
                            traps0 = grid,
                            detectpar1 = list(lambda0 = L01, sigma = sigma1),
                            detectpar2 = list(lambda0 = L02, sigma = sigma2),
                            noccasions = 1,
                            criterion = c("sum_min"), # sum_min
                            spacing_m = seq(0,5000,10))
spacing.sum.min$optimum.spacing

spacing.all.min <- optimalSpacing2G(D1 = D1, D2 = D2,
                                    traps0 = grid,
                                    detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                    detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                    noccasions = 1,
                                    criterion = c("all_min"), # all_min
                                    spacing_m = seq(0,5000,10))
spacing.all.min$optimum.spacing

#use the spacing from sum.min
grid.sum.min <- Proposed.traps(poly = GKEmsks$`GKE Polygon SF obj`, alltraps = NULL, D = NULL, sigma = NULL, 
                             lambda0 = NULL, sigma.buff = NULL, grid.spacing = 3000, 
                             criterion = 4, n.reps = 100, grid = TRUE, nT = nT)

#########################################################
#Systematic clustered designs
#first find within and between optimal spacing using strata 1 for within ad S2 for btwn
#using 2G fn so I don't need a bunch of new fns
#########################################################

grid <- make.grid(7,7,1000, detector = "count")
spacing.within <- optimalSpacing2G(D1 = D1, D2 = D1,
                                    traps0 = grid,
                                    detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                    detectpar2 = list(lambda0 = L01, sigma = sigma1),
                                    noccasions = 1,
                                    criterion = c("sum_min"), # sum_min
                                    spacing_m = seq(0,1000,200))
spacing.within$optimum.spacing
plot(spacing.within$values$crit ~ spacing.within$values$spacing, type = 'l')

grid2 <- make.grid(7,7,1000, detector = "count")
spacing.btwn <- optimalSpacing2G(D1 = D2, D2 = D2,
                                   traps0 = grid2,
                                   detectpar1 = list(lambda0 = L02, sigma = sigma2),
                                   detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                   noccasions = 1,
                                   criterion = c("sum_min"), # sum_min
                                   spacing_m = seq(0,10000,200))
spacing.btwn$optimum.spacing
plot(spacing.btwn$values$crit ~ spacing.btwn$values$spacing, type = 'l')

clusters.10 <- Proposed.traps(poly = GKEmsks$`GKE Polygon SF obj`, alltraps = NULL, D = NULL, sigma = NULL, 
                  lambda0 = NULL, sigma.buff = NULL, grid.spacing = 4000, 
                  criterion = 4, n.reps = 5, grid = TRUE, nT = 10)

# Wrapper that generates one full set of clustered trap coordinates
trap_list <- lapply(clusters.10, function(centroids) {
  make_cluster_traps(centroids, spacing = 1000, n_rows = 2, n_cols = 2)
})

trap_secr_list <- lapply(trap_list, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

Enrm(D = D1, trap_secr_list[[1]], GKEmsks$`GKE Mask`, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap_secr_list[[1]], GKEmsks$`GKE Mask`, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

####################################################

grid.designs <- list("600 m (2G opt)" = grid.sum.min, "Cluster" = trap_secr_list)
save(grid.designs, file = "GKE/GridDesigns.RData")

##########################################################
#Both strata
load("GKE/Cluster/GA4BothStrataPars.RData")



test <- Proposed.traps(poly = NULL, alltraps = GKE_trap_locs_red, mask.buff = GKE_msk, D = D1, sigma = sigma1,
                       lambda0 = L01, sigma.buff = sigma1, grid.spacing = NULL, criterion = 4, 
                       n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 50)




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

save(GA.designs, file = "GKE/GADesigns.RData")
