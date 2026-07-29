#July 2026, designs as we vary the diff btwn the strata (K)

library(secrdesign)
library(kofnGA)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

rm(list=ls())

source("functions2.R")
source('optimalSpacing2G/optimalSpacing2G.R')

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

#determine optimal spacing using 2 goups
grid <- make.grid(7,7,1000, detector = "count")
spacing.meanCV <- optimalSpacing2G(D1 = D1, D2 = D2,
                                   traps0 = grid,
                                   detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                   detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                   noccasions = 1,
                                   criterion = c("min_mean"), # mean CV
                                   spacing_m = seq(0,10000,100))
spacing.meanCV$optimum.spacing

spacing.minmax<- optimalSpacing2G(D1 = D1, D2 = D2,
                                 traps0 = grid,
                                 detectpar1 = list(lambda0 = L01, sigma = sigma1),
                                 detectpar2 = list(lambda0 = L02, sigma = sigma2),
                                 noccasions = 1,
                                 criterion = c("min_max"), # sum_min
                                 spacing_m = seq(0,10000,100))
spacing.minmax$optimum.spacing


#generate 500 grid designs for each value of K
#using the geometric mean approach
GridKEval <- function(polygon, sigma1, K, nT, nreps){
  
  GridK.list <- vector("list", K-1)
  
  for (k in 2:K){
    sigma2 <- k * sigma1
    h = sigma1^0.5 * sigma2^0.5
    h = round(h/10)*10
    GridK.list[[k-1]] <- Proposed.traps(poly = polygon, alltraps = NULL, D = NULL, sigma = NULL, 
                      lambda0 = NULL, sigma.buff = NULL, grid.spacing = h, 
                      criterion = 4, n.reps = nreps, grid = TRUE, nT = nT)
    names(GridK.list)[k - 1] <- paste0("K", k)
  }
  return(GridK.list)
}

#for 40 traps
gridK.40 <- GridKEval(traplocs.sf, sigma1 = 200, K = 14, nT = 40, nreps = 500)

save(gridK.40, file = "Cluster/Sims/GridK40.RData")

#spacings used
h.vec <- NULL ; sigma1 = 200
L0.vec <- NULL ; L01 = 2
D.vec <- NULL ; D1 = 0.05
Sig.vec <- NULL

for (k in 2:14){
  sigma2 <- k * sigma1
  h = sigma1^0.5 * sigma2^0.5
  h = round(h/10)*10
  h.vec[k-1] <- h
  
  Sig.vec[k-1] <- sigma2
  L0.vec[k-1] <- L01 / k
  D.vec[k-1] <- D1 / k
}

h.vec
Sig.vec
L0.vec
D.vec
