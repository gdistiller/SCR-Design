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
sigma1 = 200 ; sigma2 = 3000  ; sigma <- c(sigma1, sigma2)
D1 = 0.05 ; D2 = D1/15  ; D <- c(D1, D2)
L01 <- 2 ; L02 <- L01/15
nT1 <- 40 ; nT2 <- 120 ; nreps = 500

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
                                    spacing_m = seq(0,5000,100))
spacing.sum.min$optimum.spacing

spacing.all.min <- optimalSpacing2G(D1 = D1, D2 = D2,
                                    traps0 = grid,
                                    detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                    detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                    noccasions = 1,
                                    criterion = c("all_min"), # all_min
                                    spacing_m = seq(0,5000,10))
spacing.all.min$optimum.spacing

#use the spacing from all.min of 700 m
grid.sum.min.40 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                               lambda0 = NULL, sigma.buff = NULL, grid.spacing = 700, 
                               criterion = 4, n.reps = nreps, grid = TRUE, nT = nT1)

grid.sum.min.120 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                                  lambda0 = NULL, sigma.buff = NULL, grid.spacing = 700, 
                                  criterion = 4, n.reps = nreps, grid = TRUE, nT = nT2)

#800
grid.800.40 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                           lambda0 = NULL, sigma.buff = NULL, grid.spacing = 800, 
                           criterion = 4, n.reps = nreps, grid = TRUE, nT = nT1)

grid.800.120 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = 800, 
                              criterion = 4, n.reps = nreps, grid = TRUE, nT = nT2)

#1600
grid.1600.40 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                           lambda0 = NULL, sigma.buff = NULL, grid.spacing = mean(sigma), 
                           criterion = 4, n.reps = nreps, grid = TRUE, nT = nT1)

grid.1600.120 <- Proposed.traps(poly = traplocs.sf, alltraps = NULL, D = NULL, sigma = NULL, 
                               lambda0 = NULL, sigma.buff = NULL, grid.spacing = mean(sigma), 
                               criterion = 4, n.reps = nreps, grid = TRUE, nT = nT2)

#########################################################
#Systematic clustered designs
#Going with 10 clusters of 2x2
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
                                   spacing_m = seq(0,1000,100))
spacing.within$optimum.spacing

grid2 <- make.grid(7,7,1000, detector = "count")
spacing.btwn <- optimalSpacing2G(D1 = D2, D2 = D2,
                                 traps0 = grid2,
                                 detectpar1 = list(lambda0 = L02, sigma = sigma2),
                                 detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                 noccasions = 1,
                                 criterion = c("sum_min"), # sum_min
                                 spacing_m = seq(0,10000,100))
spacing.btwn$optimum.spacing

#generate a between-cluster grid for cluster centroids, so nT here is 10 or 30
#slightly increase the buffer, for 2x2 grids with 500 m spacing this is sqrt(2*250^2) or 353
res.objs1 <- create.extent(sigma = 3000, buff.factor = 3.12, res = 200)
mask1 <- res.objs1[[1]]
clust.locs1 <- res.objs1[[2]]
clustlocs.sf1 <- res.objs1[[3]]

clusters.os.40 <- Proposed.traps(poly = clustlocs.sf1, alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = spacing.btwn$optimum.spacing + spacing.within$optimum.spacing, 
                              criterion = 4, n.reps = 500, grid = TRUE, nT = 10)

clusters.os.120 <- Proposed.traps(poly = clustlocs.sf1, alltraps = NULL, D = NULL, sigma = NULL, 
                                 lambda0 = NULL, sigma.buff = NULL, grid.spacing = spacing.btwn$optimum.spacing + spacing.within$optimum.spacing, 
                                 criterion = 4, n.reps = 500, grid = TRUE, nT = 30)

plot(mask, axes = T)
plot(clusters.os.40[[1]], add = T)
plot(clusters.os.120[[1]], add = T)

# Wrapper that generates one full set of clustered trap coordinates
trap.list.os.40 <- lapply(clusters.os.40, function(centroids) {
  make_cluster_traps(centroids, spacing = 500, n_rows = 2, n_cols = 2)
})

trap.list.os.120 <- lapply(clusters.os.120, function(centroids) {
  make_cluster_traps(centroids, spacing = 500, n_rows = 2, n_cols = 2)
})

trap.secr.list.os.40 <- lapply(trap.list.os.40, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

trap.secr.list.os.120 <- lapply(trap.list.os.120, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

plot(mask1, axes = T)
plot(trap.secr.list.os.40[[10]], add = T)
plot(trap.secr.list.os.120[[10]], add = T)

Enrm(D = D1, trap.secr.list.os.40[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap.secr.list.os.40[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

Enrm(D = D1, trap.secr.list.os.120[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap.secr.list.os.120[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

#again with 2 sigma spacing (400 / 6000)
#slightly increase the buffer, for 2x2 grids with 500 m spacing this is sqrt(2*200^2) or 282.84
res.objs2 <- create.extent(sigma = 3000, buff.factor = 3.095, res = 200)
mask2 <- res.objs2[[1]]
clust.locs2 <- res.objs2[[2]]
clustlocs.sf2 <- res.objs2[[3]]

clusters.2sig.40 <- Proposed.traps(poly = clustlocs.sf2, alltraps = NULL, D = NULL, sigma = NULL, 
                                 lambda0 = NULL, sigma.buff = NULL, grid.spacing = 6400, 
                                 criterion = 4, n.reps = 500, grid = TRUE, nT = 10)

clusters.2sig.120 <- Proposed.traps(poly = clustlocs.sf2, alltraps = NULL, D = NULL, sigma = NULL, 
                                  lambda0 = NULL, sigma.buff = NULL, grid.spacing = 6400, 
                                  criterion = 4, n.reps = 500, grid = TRUE, nT = 30)

plot(mask2, axes = T)
plot(clusters.2sig.40[[1]], add = T)
plot(clusters.2sig.120[[1]], add = T)

# Wrapper that generates one full set of clustered trap coordinates
trap.list.2sig.40 <- lapply(clusters.2sig.40, function(centroids) {
  make_cluster_traps(centroids, spacing = 500, n_rows = 2, n_cols = 2)
})

trap.list.2sig.120 <- lapply(clusters.2sig.120, function(centroids) {
  make_cluster_traps(centroids, spacing = 500, n_rows = 2, n_cols = 2)
})

trap.secr.list.2sig.40 <- lapply(trap.list.2sig.40, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

trap.secr.list.2sig.120 <- lapply(trap.list.2sig.120, function(traps) {
  read.traps(data = traps[, c("x", "y")], detector = "count")
})

plot(mask2, axes = T)
plot(trap.secr.list.os.40[[10]], add = T)
plot(trap.secr.list.os.120[[10]], add = T)

Enrm(D = D1, trap.secr.list.os.40[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap.secr.list.os.40[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

Enrm(D = D1, trap.secr.list.os.120[[1]], mask, detectpar = list(lambda0 = L01, sigma = sigma1), noccasions = 1)
Enrm(D = D2, trap.secr.list.os.120[[1]], mask, detectpar = list(lambda0 = L02, sigma = sigma2), noccasions = 1)

#save designs
grid.designs.40 <- list("700 m (2G opt)" = grid.sum.min.40, "800 m" = grid.800.40, "Avg sigma" = grid.1600.40, 
                     "Cluster (opt)" = trap.secr.list.os.40, "Cluster (2 sig)" = trap.secr.list.2sig.40)

grid.designs.120 <- list("700 m (2G opt)" = grid.sum.min.120, "800 m" = grid.800.120, "Avg sigma" = grid.1600.120, 
                        "Cluster (opt)" = trap.secr.list.os.120, "Cluster (2 sig)" = trap.secr.list.2sig.120)

save(grid.designs.40, file = "15FoldScen/GridDesigns40.RData")
save(grid.designs.120, file = "15FoldScen/GridDesigns120.RData")

##########################################################
#lacework
##########################################################
# function to subset lacework designs to contain N points
source("crop_to_n.R")

# set spacings according to 2 * sigma
smallspacing <- 400
bigspacing <- 6000

res.objs <- create.extent(sigma = 3000, buff.factor = 3, res = 200)
mask <- res.objs[[1]]
trap.locs <- res.objs[[2]]
traplocs.sf <- res.objs[[3]]

region = data.frame(x = c(min(mask$x), min(mask$x), max(mask$x), max(mask$x)),
                    y = c(min(mask$y), max(mask$y), max(mask$y), min(mask$y)))

# make the design 
lwrot <- 45
lw <- make.lacework(region = region, 
                    spacing = c(bigspacing, smallspacing),  
                    rotate = lwrot,
                    radius = 1000,
                    detector = "count", 
                    keep.design = TRUE)

## Manually remove some traps so end up with the same number of traps as other designs
lw.40.1 <- crop_to_n(lw$x, lw$y, N = 40, xy_ratio = 1)
lw.40.2 <- crop_to_n(lw$x, lw$y, N = 40, xy_ratio = 2)

lw.120.1 <- crop_to_n(lw$x, lw$y, N = 120, xy_ratio = 1)
lw.120.2 <- crop_to_n(lw$x, lw$y, N = 120, xy_ratio = 2)

lw.40.1 <- read.traps(data = lw.40.1$points, detector = "count")
lw.40.2 <- read.traps(data = lw.40.2$points, detector = "count")

lw.120.1 <- read.traps(data = lw.120.1$points, detector = "count")
lw.120.2 <- read.traps(data = lw.120.2$points, detector = "count")

plot(mask)
plot(lw.40.1, add = TRUE)

plot(mask)
plot(lw.40.2, add = TRUE)

plot(mask)
plot(lw.120.1, add = TRUE)

plot(mask)
plot(lw.120.2, add = TRUE)

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