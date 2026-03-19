#28 March 2025, run again in November to get 500 reps
#March 26, adding optimal spacing, and clustered systematic
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

#create grid designs using optimal spacing, sum_min,  all_min critertia for two classes
#result does not seem sensitive to spacing of base grid
grid <- make.grid(7,7,1000, detector = "count")
spacing.sum.min <- optimalSpacing2G(D1 = D1, D2 = D2,
                            traps0 = grid,
                            detectpar1 = list(lambda0 = L01, sigma = sigma1),
                            detectpar2 = list(lambda0 = L02, sigma = sigma2),
                            noccasions = 1,
                            criterion = c("sum_min"), # sum_min
                            spacing_m = seq(100,5000,100))
spacing.sum.min$optimum.spacing

spacing.all.min <- optimalSpacing2G(D1 = D1, D2 = D2,
                                    traps0 = grid,
                                    detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                    detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                    noccasions = 1,
                                    criterion = c("all_min"), # all_min
                                    spacing_m = seq(100,5000,100))
spacing.all.min$optimum.spacing

#use the spacing from all.min of 700 m
grid.sum.min <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                             lambda0 = NULL, sigma.buff = NULL, grid.spacing = 700, 
                             criterion = 4, n.reps = nreps, grid = TRUE, nT = nT)

#800
grid.800 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                           lambda0 = NULL, sigma.buff = NULL, grid.spacing = 800, 
                           criterion = 4, n.reps = nreps, grid = TRUE, nT = nT)

#1600
grid.1600<- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = mean(sigma), 
                              criterion = 4, n.reps = nreps, grid = TRUE, nT = nT)

grid.designs <- list("700 m (opt)" = grid.sum.min, "800 m" = grid.800, "Avg sigma" = grid.1600)

save(grid.designs, file = "15FoldScen/GridDesigns.RData")
save(res.objs, file = "15FoldScen/SCRObjs.RData")

#########################################################
#Systematic clustered designs
#first find within and between optimal spacing using strata 1 for within a S2 for btwn
#using 2G fn so I dont need a bunch of new fns
#also some issue with dfcast (similar to secr_valid.detecfn)
grid <- make.grid(7,7,100, detector = "count")
spacing.within <- optimalSpacing2G(D1 = D1, D2 = D1,
                                    traps0 = grid,
                                    detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                    detectpar2 = list(lambda0 = L01, sigma = sigma1),
                                    noccasions = 1,
                                    criterion = c("sum_min"), # sum_min
                                    spacing_m = seq(200,1000,200))
spacing.within$optimum.spacing

grid2 <- make.grid(7,7,1000, detector = "count")
spacing.btwn <- optimalSpacing2G(D1 = D2, D2 = D2,
                                   traps0 = grid2,
                                   detectpar1 = list(lambda0 = L02, sigma = sigma2),
                                   detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                   noccasions = 1,
                                   criterion = c("sum_min"), # sum_min
                                   spacing_m = seq(1000,10000,200))
spacing.btwn$optimum.spacing

#generate a between-cluster grid for cluster centroids
#slightly increase the buffer, for 2x2 grids this is sqrt(2*300^2) or 424.26
res.objs <- create.extent(sigma = 3000, buff.factor = 3.15, res = 200)
mask <- res.objs[[1]]
clust.locs <- res.objs[[2]]
clustlocs.sf <- res.objs[[3]]

clusters.10 <- Proposed.traps(poly = clustlocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                  lambda0 = NULL, sigma.buff = NULL, grid.spacing = spacing.btwn$optimum.spacing + spacing.within$optimum.spacing, 
                  criterion = 4, n.reps = 500, grid = TRUE, nT = 10)

plot(mask, axes = T)
plot(clusters.10[[1]], add = T)

# Wrapper that generates one full set of clustered trap coordinates
trap_list <- lapply(clusters.10, function(centroids) {
  make_cluster_traps(centroids, spacing = 600, n_rows = 2, n_cols = 2)
})

trap_secr_list <- lapply(trap_list, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

plot(mask, axes = T)
plot(trap_secr_list[[10]], add = T)

Enrm(D = D1, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

#again with 400 / 1500 spacing
clusters.10 <- Proposed.traps(poly = clustlocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = 1900, 
                              criterion = 4, n.reps = 2, grid = TRUE, nT = 10)

# Wrapper that generates one full set of clustered trap coordinates
trap_list <- lapply(clusters.10, function(centroids) {
  make_cluster_traps(centroids, spacing = 400, n_rows = 2, n_cols = 2)
})

trap_secr_list <- lapply(trap_list, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

Enrm(D = D1, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

#again with 400 / 6000 spacing
clusters.10 <- Proposed.traps(poly = clustlocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = 6400, 
                              criterion = 4, n.reps = 2, grid = TRUE, nT = 10)

# Wrapper that generates one full set of clustered trap coordinates
trap_list <- lapply(clusters.10, function(centroids) {
  make_cluster_traps(centroids, spacing = 400, n_rows = 2, n_cols = 2)
})

trap_secr_list <- lapply(trap_list, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

Enrm(D = D1, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

#again with 400 / 3000 spacing
clusters.10 <- Proposed.traps(poly = clustlocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = 3400, 
                              criterion = 4, n.reps = 2, grid = TRUE, nT = 10)

# Wrapper that generates one full set of clustered trap coordinates
trap_list <- lapply(clusters.10, function(centroids) {
  make_cluster_traps(centroids, spacing = 400, n_rows = 2, n_cols = 2)
})

trap_secr_list <- lapply(trap_list, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

Enrm(D = D1, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap_secr_list[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

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
