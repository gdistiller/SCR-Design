#Enrm and data summaries for large reg area

library(ggplot2)
library(secr)
library(secrdesign)
library(kofnGA)
library(dplyr)
library(tidyr)
library(sf)
library(terra)

#create mask and potential traplocs
mask.r <- rast(ncol = 120, nrow = 120, xmin = 0, xmax = 150000, ymin = 0, ymax = 150000)  #ncol/row determine the resolution
traplocs.r <- rast(ncol = 60, nrow = 60, xmin = 30000, xmax = 120000, ymin = 30000, ymax = 120000) #res of 1.5K

values(mask.r) <- 1:ncell(mask.r) #numbers each cell, seems needed
values(traplocs.r) <- 1:ncell(traplocs.r) 

#create polyon
msk.p1 <- rbind(c(0, 0), c(150000, 0), c(150000, 150000), c(0, 150000), c(0,0))
mask.v <- vect(msk.p1, type = "polygons")

traplocs.p1 <- rbind(c(30000, 30000), c(120000, 30000), c(120000, 120000), c(30000, 120000), c(30000,30000))
traplocs.v <- vect(traplocs.p1, type = "polygons")

#create mask and traps obj
#get coordinates
msk.coords <- crds(mask.r, df = T)
traploc.coords <- crds(traplocs.r, df=T)

mask <- read.mask(data = msk.coords, spacing = 1250)
trap.locs <- read.traps(data = traploc.coords, detector = "count", spacing = 1500)

# 40 traps
D4 = 0.0004 ; L0 = 3 ; sig = 3000

nT <- 40

######
#Grid, spacing of 2 km (2/3 sigma)
######

#for grids just need to pass the polygon (for trap locs)
locs <- st_as_sf(traplocs.v)
msk <- st_as_sf(mask.v)

grid.2 <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                          lambda0 = L0, sigma.buff = NULL, grid.spacing = 2000, 
                          criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

plot(mask.v)
plot(grid.2[[1]], add = T)

######
#GA
######
GA.F4 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 40)
plot(mask.v)
plot(GA.F4[[1]], add = T)

GA.F5 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask,  D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 5, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 40)
plot(mask.v)
plot(GA.F5[[1]], add = T)

#Enrm numbers
scenF <- make.scenarios (trapsindex = 1:3, detectfn = 'HHN', D = D4, lambda0 = L0, 
                         sigma = sig, noccasions = 1)

traplist <- list(Grid.40 = grid.2[[1]], GA.fill = GA.F4[[1]], GA.cluster = GA.F5[[1]])

scenarioSummary(scenF, traps = traplist, mask = mask)

test.data <- run.scenarios(20, scenF, traplist, mask, fit = FALSE, extract = identity, 
                           det.args = list(savepopn = TRUE), fit.args = list(detectfn = 'HHN'))
summary(test.data)

test.fit <- run.scenarios(20, scenF, traplist, mask, fit = TRUE, fit.args = list(detectfn = 'HHN'), ncores = 4)
