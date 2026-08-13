#Aug 2026, 2 stage approach design
#redoing with En2 for stage 2
#will be added to the design objects afterwards

library(secrdesign)
library(kofnGA)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

rm(list=ls())

source("functions2.R")
#note Ian's GAoptim fns added to functions2

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

# Use all trap locations as the sf sample frame for a spatially balanced sample

nfixed1 <- nT1/2 ; nfixed2 <- nT2/2
alltrps_sf <- sf::st_as_sf(data.frame(id = seq_len(nrow(trap.locs)),
                                      x = trap.locs$x,
                                      y = trap.locs$y),
                           coords = c("x", "y"))

alltrps_sf_poly <- sf::st_convex_hull(sf::st_union(alltrps_sf))
alltrps_poly_sf <- sf::st_sf(id = 1, geometry = alltrps_sf_poly)

#first for 40 traps
grts_out <- spsurvey::grts(alltrps_poly_sf, n_base = nfixed1, projcrs_check = FALSE)
fixed_idx <- sf::st_nearest_feature(grts_out$sites_base, alltrps_sf)
fixedtrps_sf <- trap.locs[fixed_idx, ]
remtrps_sf <- trap.locs[-fixed_idx, ]

plot(mask)
plot(trap.locs, detpar = list(pch = 1, cex = 0.4, fg = "lightgreen"), add = TRUE)
plot(fixedtrps_sf, detpar = list(fg = "black"), add = TRUE)

opt.40 <- GAoptim(mask, alltraps = remtrps_sf, fixedtraps = fixedtrps_sf, ntraps = nT1,
                detectpar = list(lambda0 = L01, sigma = sigma1),
                detectfn = 'HHN', D = D1, noccasions = 1, ngen = 10, verbose = 1)

plot(mask)
plot(trap.locs, detpar = list(pch = 1, cex = 0.4, fg = "lightgreen"), add = TRUE)
plot(opt.40$fixedtraps, detpar = list(fg = "forestgreen"), add = TRUE)
plot(opt.40$selectedtraps, add = TRUE)
plot(opt.40$optimaltraps, detpar = list(pch = 21, cex = 1.2), add = TRUE)
minnrRSE(opt.40, distribution = 'binomial')

#and 120 traps
grts_out <- spsurvey::grts(alltrps_poly_sf, n_base = nfixed2, projcrs_check = FALSE)
fixed_idx <- sf::st_nearest_feature(grts_out$sites_base, alltrps_sf)
fixedtrps_sf <- trap.locs[fixed_idx, ]
remtrps_sf <- trap.locs[-fixed_idx, ]

plot(mask)
plot(trap.locs, detpar = list(pch = 1, cex = 0.4, fg = "lightgreen"), add = TRUE)
plot(fixedtrps_sf, detpar = list(fg = "black"), add = TRUE)

opt.120 <- GAoptim(mask, alltraps = remtrps_sf, fixedtraps = fixedtrps_sf, ntraps = nT2,
                  detectpar = list(lambda0 = L01, sigma = sigma1),
                  detectfn = 'HHN', D = D1, noccasions = 1, ngen = 10, verbose = 1)

plot(mask)
plot(trap.locs, detpar = list(pch = 1, cex = 0.4, fg = "lightgreen"), add = TRUE)
plot(opt.120$fixedtraps, detpar = list(fg = "forestgreen"), add = TRUE)
plot(opt.120$selectedtraps, add = TRUE)
plot(opt.120$optimaltraps, detpar = list(pch = 21, cex = 1.2), add = TRUE)
minnrRSE(opt.120, distribution = 'binomial')

#now generate nreps designs using a function
TwoStage.40 <- Two.stage.design(msk = mask, alltraps = trap.locs, L0 = L01, Sig = sigma1, Dens = D1, ntrps = nT1, nfixed = nT1/2, ngens = 3, nreps = 1)

TwoStage.120 <- Two.stage.design(msk = mask, alltraps = trap.locs, L0 = L01, Sig = sigma1, Dens = D1, ntrps = nT2, nfixed = nT2/2, ngens = 200, nreps = 1)

#500 reps run on cluster
#I used .combine = 'c' so there are 3 elements for each rep, wil extract by element name
#1st block does stage 2 with min(n,r), redone afterwards using En2
setwd("~/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/40Traps")
setwd("C:/Users/Greg/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/40Traps")

#40 traps

load("GA2Stage40b.RData")

GA2Stage40 <- TwoStage.40.list[names(TwoStage.40.list) == "Proposed design"]

#reformat so the obj contains a list of 500 (rather than a list of lists)
GA2Stage40 <- lapply(GA2Stage40, `[[`, 1)

GA2StageDesigns <- list("40 traps" = GA2Stage40)
#created without 120 traps while that is running


#120 traps, first doing with set from ngen = 300
setwd("C:/Users/Greg/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/120Traps")

load("GA2Stage120c.RData")
GA2Stage120 <- TwoStage.120.list[names(TwoStage.120.list) == "Proposed design"]
GA2Stage120 <- lapply(GA2Stage120, `[[`, 1)

GA2StageDesigns <- list("40 traps" = GA2Stage40, "120 traps" = GA2Stage120)

save(GA2StageDesigns, file = "Cluster/Sims/GA2StageDesigns.RData")

#and again using version b with ngen = 500
setwd("C:/Users/Greg/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/120Traps")

load("GA2Stage120b.RData")
GA2Stage120 <- TwoStage.120.list[names(TwoStage.120.list) == "Proposed design"]
GA2Stage120 <- lapply(GA2Stage120, `[[`, 1)

GA2StageDesigns <- list("40 traps" = GA2Stage40, "120 traps" = GA2Stage120)

save(GA2StageDesigns, file = "Cluster/Sims/GA2StageDesignsb.RData")

#################################################################
#now with En2
setwd("~/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/40Traps")
setwd("C:/Users/Greg/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/40Traps")

#40 traps (version c for 40 and d for 120)

load("GA2Stage40c.RData")

GA2Stage40 <- TwoStage.40.list[names(TwoStage.40.list) == "Proposed design"]

#reformat so the obj contains a list of 500 (rather than a list of lists)
GA2Stage40 <- lapply(GA2Stage40, `[[`, 1)

setwd("C:/Users/Greg/OneDrive - University of Cape Town/Documents/Git/SCRDesign/15FoldScen/Cluster/ProposedDesigns/120Traps")

load("GA2Stage120d.RData")
GA2Stage120 <- TwoStage.120.list[names(TwoStage.120.list) == "Proposed design"]
GA2Stage120 <- lapply(GA2Stage120, `[[`, 1)

GA2StageDesigns <- list("40 traps" = GA2Stage40, "120 traps" = GA2Stage120)

save(GA2StageDesigns, file = "Cluster/Sims/GA2StageDesignsEn2.RData")

