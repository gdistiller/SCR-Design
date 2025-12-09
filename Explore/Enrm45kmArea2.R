#April 2025, redoing after increasing ngen
#create regular area using approx 15 * M  sigma
#core unbuffered area 27 km or about 9 sigma
#creating area as both polygon and raster object
#polygon required for randomised grid designs, raster -> traps / mask required for GA designs

library(secrdesign)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

rm(list=ls())

######################################
#function to create a mask and trap locs objects, includes SCR mask and traps obj and a SF polygon for grids
#sig.factor determines the extent ito sigma, buffer is equal to 3 sigma
#res is the resolution of both
create.extent <- function(sigma, sig.factor = 15, res = 200){
  Extent <- sig.factor * sigma ; ncolrow1 = Extent/res 
  BuffExtent1 <- 3 * sigma ; BuffExtent2 <- Extent-BuffExtent1 ; ncolrow2 = (BuffExtent2-BuffExtent1)/res
  
  mask.r <- rast(ncol = ncolrow1, nrow = ncolrow1, xmin = 0, xmax = Extent, ymin = 0, ymax = Extent)
  traplocs.r <- rast(ncol = ncolrow2, nrow = ncolrow2, xmin = BuffExtent1, xmax = BuffExtent2, ymin = BuffExtent1, ymax = BuffExtent2) 
  
  values(mask.r) <- 1:ncell(mask.r) #numbers each cell, seems needed
  values(traplocs.r) <- 1:ncell(traplocs.r) 
  
  #create polygon
  msk.p1 <- rbind(c(0, 0), c(Extent, 0), c(Extent, Extent), c(0, Extent), c(0,0))
  mask.v <- vect(msk.p1, type = "polygons")
  
  traplocs.p1 <- rbind(c(BuffExtent1, BuffExtent1), c(BuffExtent2, BuffExtent1), c(BuffExtent2, BuffExtent2), c(BuffExtent1, BuffExtent2), c(BuffExtent1,BuffExtent1))
  traplocs.v <- vect(traplocs.p1, type = "polygons")
  traplocs.sf <- sf::st_as_sf(traplocs.v)
  
  #create mask and traps obj
  #get coordinates
  msk.coords <- crds(mask.r, df = T)
  traploc.coords <- crds(traplocs.r, df=T)
  
  mask <- read.mask(data = msk.coords, spacing = res)
  trap.locs <- read.traps(data = traploc.coords, detector = "count", spacing = res)
  objs <- list("SCR Mask" = mask, "SCR trap locs" = trap.locs, "SF traps polygon" = traplocs.sf)
}

#function to get summaries for a given trap obj, D1, and L0
getEnrm <- function(trap.obj, mask.obj, Dens, L0, sigma = 200){
  scen <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = Dens, lambda0 = L0, 
                          sigma = sigma, noccasions = 1)
  scenarioSummary(scen, traps = trap.obj, mask = mask.obj)
}

######################################

#######################################
#Set values for both strata
#use constant lambda0 of 1, and 40 traps
########################################
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
D1 = 0.1 ; D2 = D1 / DiffFactor  ; D <- c(D1, D2)
lambda0 = 1
nT <- 40

################################
#45km area
################################

#create mask and potential traplocs
#using 15 * larger sigma for extent and 3 sigma for buffer
res.objs <- create.extent(sigma = 3000)
mask <- res.objs[[1]]

#grid spacing of 300 m for 40 and 100 traps
grid.300 <- Proposed.traps(poly = res.objs[[3]], alltraps = NULL, D = NULL, sigma = NULL, 
                            lambda0 = NULL, sigma.buff = NULL, grid.spacing = 300, 
                            criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.300.L <- Proposed.traps(poly = res.objs[[3]], alltraps = NULL, D = NULL, sigma = NULL, 
                           lambda0 = NULL, sigma.buff = NULL, grid.spacing = 300, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 100)

#grid spacing of 1600 m for 40 and 100 traps
grid.1600 <- Proposed.traps(poly = res.objs[[3]], alltraps = NULL, D = NULL, sigma = NULL, 
                            lambda0 = NULL, sigma.buff = NULL, grid.spacing = 1600, 
                            criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.1600.L <- Proposed.traps(poly = res.objs[[3]], alltraps = NULL, D = NULL, sigma = NULL, 
                              lambda0 = NULL, sigma.buff = NULL, grid.spacing = 1600, 
                              criterion = 4, n.reps = 1, grid = TRUE, nT = 100)

#get Enrm summaries
#spacing of 300

getEnrm(grid.300, mask, Dens = 0.01, L0 = 1.5,sigma = 200)
getEnrm(grid.300, mask, Dens = 0.01/15, L0 = 1.5,sigma = 3000)

getEnrm(grid.300, mask, Dens = 0.05, L0 = 1.5,sigma = 200)
getEnrm(grid.300, mask, Dens = 0.05/15, L0 = 1.5,sigma = 3000)

#spacing of 1600

getEnrm(grid.1600, mask, Dens = 0.01, L0 = 1.5,sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.01/15, L0 = 1.5,sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.03, L0 = 1.5,sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.03/15, L0 = 1.5,sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.05, L0 = 1.5,sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.05/15, L0 = 1.5,sigma = 3000)

getEnrm(grid.1600.L, mask, Dens = 0.1, L0 = 1.5,sigma = 200)
getEnrm(grid.1600.L, mask, Dens = 0.1/15, L0 = 1.5,sigma = 3000)

#and now fixing L0 = 1
getEnrm(grid.300, mask, Dens = 0.01, L0 = 1, sigma = 200)
getEnrm(grid.300, mask, Dens = 0.01/15, L0 = 1, sigma = 3000)

getEnrm(grid.300, mask, Dens = 0.03, L0 = 1,sigma = 200)
getEnrm(grid.300, mask, Dens = 0.03/15, L0 = 1,sigma = 3000)

getEnrm(grid.300, mask, Dens = 0.05, L0 = 1,sigma = 200)
getEnrm(grid.300, mask, Dens = 0.05/15, L0 = 1,sigma = 3000)

getEnrm(grid.300.L, mask, Dens = 0.01, L0 = 1, sigma = 200)
getEnrm(grid.300.L, mask, Dens = 0.01/15, L0 = 1, sigma = 3000)

getEnrm(grid.300.L, mask, Dens = 0.03, L0 = 1,sigma = 200)
getEnrm(grid.300.L, mask, Dens = 0.03/15, L0 = 1,sigma = 3000)

getEnrm(grid.300.L, mask, Dens = 0.05, L0 = 1,sigma = 200)
getEnrm(grid.300.L, mask, Dens = 0.05/15, L0 = 1,sigma = 3000)

#spacing of 1600
getEnrm(grid.1600, mask, Dens = 0.01, L0 = 1, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.01/15, L0 = 1, sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.03, L0 = 1,sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.03/15, L0 = 1,sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.05, L0 = 1,sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.05/15, L0 = 1,sigma = 3000)

getEnrm(grid.1600.L, mask, Dens = 0.01, L0 = 1, sigma = 200)
getEnrm(grid.1600.L, mask, Dens = 0.01/15, L0 = 1, sigma = 3000)

getEnrm(grid.1600.L, mask, Dens = 0.03, L0 = 1,sigma = 200)
getEnrm(grid.1600.L, mask, Dens = 0.03/15, L0 = 1,sigma = 3000)

getEnrm(grid.1600.L, mask, Dens = 0.05, L0 = 1,sigma = 200)
getEnrm(grid.1600.L, mask, Dens = 0.05/15, L0 = 1,sigma = 3000)

#and now varying L0 as well, with higher D and spacing of 300 and 40 traps
L01 <- 1.5 ; L02 <- L01/15

getEnrm(grid.300, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.300, mask, Dens = 0.05/15, L0 = L02,sigma = 3000)

#spacing = 1600
getEnrm(grid.1600, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.1, L0 = L01, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.1/15, L0 = L02, sigma = 3000)

#varying L0 with larger diff, with 1600 spacing
L01 <- 2 ; L02 <- L01/20

getEnrm(grid.1600, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.1, L0 = L01, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.1/15, L0 = L02, sigma = 3000)

#varying L0 with higher L01 and 15 fold diff
L01 <- 2 ; L02 <- L01/15

getEnrm(grid.1600, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

getEnrm(grid.1600, mask, Dens = 0.1, L0 = L01, sigma = 200)
getEnrm(grid.1600, mask, Dens = 0.1/15, L0 = L02, sigma = 3000)

#trying for spacing of 800m
DiffFactor <- 20
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor
nT <- 40

#grid spacing of 300 m for 40 and 100 traps
grid.800 <- Proposed.traps(poly = res.objs[[3]], alltraps = NULL, D = NULL, sigma = NULL, 
                           lambda0 = NULL, sigma.buff = NULL, grid.spacing = 800, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

getEnrm(grid.800, mask, Dens = 0.01, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.01/15, L0 = L02, sigma = 3000)

getEnrm(grid.800, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

getEnrm(grid.800, mask, Dens = 0.1, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.1/15, L0 = L02, sigma = 3000)

getEnrm(grid.800, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.05/20, L0 = L02, sigma = 3000)

getEnrm(grid.800, mask, Dens = 0.1, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.1/20, L0 = L02, sigma = 3000)

#And again with 15 fold diff in L0
L01 <- 2 ; L02 <- L01/15

getEnrm(grid.800, mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

getEnrm(grid.800, mask, Dens = 0.1, L0 = L01, sigma = 200)
getEnrm(grid.800, mask, Dens = 0.1/15, L0 = L02, sigma = 3000)

##############################################################
#for GA designs, proposed designs run on cluster
#redone after increasing ngen

#Set values for both strata, 165 fold diff with L01 = 2
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor ; L0 <- c(L01, L02)
D1 <- 0.05 ; D2 <- D1 / DiffFactor ; D <- c(D1, D2)
nT <- 40

#plot realisations of the GA designs and get Enrm

setwd("C:/Users/Greg/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen")

#get mask
load("SCRObjs.RData")
mask <- res.objs[[1]]

#load GA4/5 with strata 1
load("Cluster/GAStrata1Pars2.RData")

par(mfrow=c(1,2))
plot(mask, axes = T)
plot(GA4.Strata1Pars2.40.list[[1]], add = T)

getEnrm(GA4.Strata1Pars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA4.Strata1Pars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

plot(mask, axes = T)
plot(GA5.Strata1Pars2.40.list[[1]], add = T)

getEnrm(GA5.Strata1Pars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA5.Strata1Pars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

#load GA4/5 with strata 2
load("Cluster/GAStrata2Pars2.RData")

par(mfrow=c(1,2))
plot(mask, axes = T)
plot(GA4.Strata2Pars2.40.list[[1]], add = T)

getEnrm(GA4.Strata2Pars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA4.Strata2Pars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

plot(mask, axes = T)
plot(GA5.Strata2Pars2.40.list[[1]], add = T)

getEnrm(GA5.Strata2Pars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA5.Strata2Pars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

#load GA4/5 with avg strata  
load("Cluster/GAAvgStrataPars2.RData")

par(mfrow=c(1,2))
plot(mask, axes = T)
plot(GA4.AvgStrataPars2.40.list[[1]], add = T)

getEnrm(GA4.AvgStrataPars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA4.AvgStrataPars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

plot(mask, axes = T)
plot(GA5.AvgStrataPars2.40.list[[1]], add = T)

getEnrm(GA5.AvgStrataPars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA5.AvgStrataPars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

#load GA4/5 with both strata  
load("Cluster/GABothStrataPars2a.RData")
load("Cluster/GABothStrataPars2b.RData")

par(mfrow=c(1,2))
plot(mask, axes = T)
plot(GA4.BothStrataPars2.40.list[[1]], add = T)

getEnrm(GA4.BothStrataPars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA4.BothStrataPars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

plot(mask, axes = T)
plot(GA5.BothStrataPars2.40.list[[1]], add = T)

getEnrm(GA5.BothStrataPars2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA5.BothStrataPars2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

#load GA4/5 with both strata  (max = T)
load("Cluster/GABothStrataParsMax2a.RData")
load("Cluster/GABothStrataParsMax2b.RData")

par(mfrow=c(1,2))
plot(mask, axes = T)
plot(GA4.BothStrataParsMax2.40.list[[1]], add = T)

getEnrm(GA4.BothStrataParsMax2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA4.BothStrataParsMax2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

plot(mask, axes = T)
plot(GA5.BothStrataParsMax2.40.list[[1]], add = T)

getEnrm(GA5.BothStrataParsMax2.40.list[[1]], mask, Dens = 0.05, L0 = L01, sigma = 200)
getEnrm(GA5.BothStrataParsMax2.40.list[[1]], mask, Dens = 0.05/15, L0 = L02, sigma = 3000)

#########################################################
