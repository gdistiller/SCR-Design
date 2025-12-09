#Sept 2023, misc code for SCR design project

#################################################################################
#settling on base setup for realistic leopard data
#using 1.5 sigma for trap spacing
library(secr) ; library(secrdesign)

#SCR parameters
MsigmaFactor = 2 #male sigma double female sigma
MDFactor = 0.25  #male D 1/4 female D

sigmaF = 3000 ; sigmaM = MsigmaFactor * sigmaF
DF = 0.00015 ; DM = MDFactor * DF
lambda0F = 5 ; lambda0M = 4.5

trapsF <- make.grid(7,7, detector = 'count', spacing = 1.5*sigmaF)
trapsM <- make.grid(7,7, detector = 'count', spacing = 1.5*sigmaM)

scen.M <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = DM,
                          lambda0 = lambda0M, sigma = sigmaM, noccasions = 1)
scen.F <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = DF,
                          lambda0 = lambda0F, sigma = sigmaF, noccasions = 1)
scenarioSummary(scen.F, trapsF)
scenarioSummary(scen.M, trapsM)

#################################################################################
#extract sim results

estimateSummary(test, 'D',  format = 'data.frame', cols = c(1,2,6:8))
estimateSummary(test, 'lambda0',  format = 'data.frame',cols = c(1,2,6:8))
estimateSummary(test, 'sigma',  format = 'data.frame',cols = c(1,2,6:8))

#for new function that returns data summaries, need to go:
estimateSummary(sim.results1.M[[3]], 'D',  format = 'data.frame', cols = c(1,2,6:8))
estimateSummary(sim.results1.M[[3]], 'lambda0',  format = 'data.frame', cols = c(1,2,6:8))
estimateSummary(sim.results1.M[[3]], 'sigma',  format = 'data.frame', cols = c(1,2,6:8))

#################################################################################
#####bits of code
red_factor <- cellsize[1] / attr(msk, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")

#mask and msk both the same thing? mask objects
mask <- secr::raster(msk)
mask_df <- data.frame(coordinates(mask))
mask <- read.mask(data = mask_df)
plot(mask, dots=F)

plot(msk)
plot(opt$optimaltraps, add = TRUE)
minnrRSE(opt, distribution = 'binomial')

buffer_mult <- ceiling(bufferF / cellsizeF)
newmask_df <- expand.grid(x = seq(from = min(mask_df$x) - buffer_mult * cellsizeF, 
                                  to = max(mask_df$x) + buffer_mult * cellsizeF, 
                                  by = cellsizeF),
                          y = seq(from = min(mask_df$y) - buffer_mult * cellsizeF, 
                                  to = max(mask_df$y) + buffer_mult * cellsizeF, 
                                  by = cellsizeF)) 
newmask_df <- newmask_df %>% 
  mutate(keep = (apply(e2dist(as.matrix(newmask_df), as.matrix(mask_df)), 1, min) <= bufferF)) %>%
  filter(keep == TRUE) %>% dplyr::select(-keep)
newmask <- read.mask(data = newmask_df)
plot(newmask,axes=T)

# create grid of possible trap locations
cellsize <- 2/3 * sigma
red_factor <- cellsize[1] / attr(msk, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")

alltraps_df <- data.frame(coordinates(data.frame(alltrps)))
alltraps <- as.matrix(alltraps_df)[,c(1,2)]
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)  


mask.raster <- secr::raster(mask)
mask_df <- data.frame(coordinates(mask.raster))
mask <- read.mask(data = mask_df)

# create grid of possible trap locations
alltraps_df <- data.frame(coordinates(data.frame(alltraps)))
alltraps <- as.matrix(alltraps_df)

###################################
#testing missing designs
MsigmaFactor = 2.5 #male sigma 2.5 female sigma
MDFactor = 0.25  #male D 1/4 female D

sigmaF = 2000 ; sigmaM = MsigmaFactor * sigmaF
DF = 0.00015 ; DM = MDFactor * DF 
lambda0F = 5 ; lambda0M = 4.5 

alltrps <- make.grid(nx = 10, ny = 10,  spacing = 2000, detector = "count")
mask <- make.mask(alltrps, buffer = 0)
mask.buff <- make.mask(alltrps, buffer = sigmaF)

n.reps = 1

set.seed(1234)
grid.traps <- Proposed.traps(mask = NULL, alltraps = alltrps, nT = 8, D = NULL, sigma = sigmaF, lambda0M,  n.reps = n.reps, grid = TRUE)

plot(mask.buff, axes = T)

for (i in 1:n.reps){
  plot(mask.buff, axes = T)
  plot(GA.traps[[i]],add=T)
  Sys.sleep(0.5)
}

#GA
GA.traps <- Proposed.traps(mask = NULL, alltraps = alltrps, nT = 8, D = NULL, sigma = sigmaF, lambda0M,  n.reps = n.reps, grid = FALSE)

plot(mask.buff, axes = T)
plot(GA.traps[[1]],add=T)

for (i in 1:n.reps){
  plot(mask.buff, axes = T)
  plot(GA.traps[[i]],add=T)
  Sys.sleep(0.5)
}

###############################################
#get Enrm for grid that matches Sims1
#SCR parameters that match 
MsigmaFactor = 2 #male sigma 2.5 female sigma
MDFactor = 0.25  #male D 1/4 female D

sigmaF = 3000 ; sigmaM = MsigmaFactor * sigmaF
DF = 0.00015 ; DM = MDFactor * DF 
lambda0F = 5 ; lambda0M = 4.5 

nT <- 49

alltrpsF <- make.grid(nx = 50, ny = 50, origin = c(1000,1000), spacing = 2/3 *sigmaF, detector = "count")
alltrpsM <- make.grid(nx = 50, ny = 50, origin = c(1000,1000), spacing = 2/3 *sigmaM, detector = "count")
alltrpsBoth <- make.grid(nx = 50, ny = 50, origin = c(1000,1000), spacing = mean(c(2/3 *sigmaF, 2/3 * sigmaM)), detector = "count")

#get Enrm for GA designs
GA.F = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DF,DM), sigma = c(sigmaF,sigmaM), 
                  lambda0 = c(lambda0F, lambda0M), nreps = 1, grid = FALSE)
GA.M = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DM,DF), sigma = c(sigmaM,sigmaF), lambda0 = c(lambda0M, lambda0F), 
                  nreps = 1, grid = FALSE)

summary(GA.F[,9:12])
summary(GA.M[,9:12])

#get Enrm for grid designs
grids.F = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DF,DM), sigma = c(sigmaF,sigmaM), lambda0 = lambda0F, nreps = 5, grid = TRUE)
grids.M = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DM,DF), sigma = c(sigmaM,sigmaF), lambda0 = lambda0M, nreps = 5, grid = TRUE)
summary(grids.F[,9:12])
summary(grids.M[,9:12])


#CLT type values
nT = 80 ; DF = 0.0001 ; DM= 0.000075 ; sigmaF = 2000 ; sigmaM = 5000 ; lambda0F = 3 ; lambda0M = 3
grids.F = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DF,DM), sigma = c(sigmaF,sigmaM), lambda0 = lambda0F, nreps = 3, grid = TRUE)
grids.M = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DM,DF), sigma = c(sigmaM,sigmaF), lambda0 = lambda0M, nreps = 3, grid = TRUE)
summary(grids.F[,9:12])
summary(grids.M[,9:12])

#increase D
nT = 80 ; DF = 0.0002 ; DM= 0.0001 ; sigmaF = 2000 ; sigmaM = 5000 ; lambda0F = 3 ; lambda0M = 3
grids.F = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DF,DM), sigma = c(sigmaF,sigmaM), lambda0 = lambda0F, nreps = 5, grid = TRUE)
grids.M = check.Enrm(alltraps = alltrpsF, nT = nT, D = c(DM,DF), sigma = c(sigmaM,sigmaF), lambda0 = lambda0M, nreps = 5, grid = TRUE)
summary(grids.F[,9:12])
summary(grids.M[,9:12])


alltrpsFlarge <- make.grid(nx = 100, ny = 100, origin = c(1000,1000), spacing = 2/3 *sigmaF, detector = "count")
grids.F = check.Enrm(alltraps = alltrpsFlarge, nT = nT, D = c(DF,DM), 
                     sigma = c(sigmaF,sigmaM), lambda0 = 3, nreps = 1, grid = TRUE)



grids.F = check.Enrm(alltraps = alltrpsF, nT = 90, D = c(0.00015,DM), sigma = c(sigmaF,sigmaM), lambda0 = 3, nreps =3, grid = TRUE)
grids.M = check.Enrm(alltraps = alltrpsF, nT = 90, D = c(DM,DF), sigma = c(sigmaM,sigmaF), lambda0 = 3, nreps = 3, grid = TRUE)

summary(grids.F[,9:12])
summary(grids.M[,9:12])

grids.F = check.Enrm(alltraps = alltrpsF, nT = 90, D = c(0.0002,0.0001), sigma = c(sigmaF,sigmaM), lambda0 = 2, nreps = 3, grid = TRUE)
grids.M = check.Enrm(alltraps = alltrpsF, nT = 90, D = c(0.0001, 0.00015), sigma = c(sigmaM,sigmaF), lambda0 = 1.5, nreps = 3, grid = TRUE)

summary(grids.F[,9:12])
summary(grids.M[,9:12])

grids.F = check.Enrm(alltraps = alltrpsF, nT = 90, D = c(0.0003,DM), sigma = c(sigmaF,sigmaM), lambda0 = 1.5, nreps = 1, grid = TRUE)

grids.M = check.Enrm(alltraps = alltrpsM, nT = 90, D = c(0.000125, 0.0003), sigma = c(sigmaM,sigmaF), lambda0 = 1.5, nreps = 1, grid = TRUE)



#first getting proposed designs, then Enr
alltrpsF <- make.grid(nx = 50, ny = 50, spacing = 2/3 *sigmaF, detector = "count")
mask <- make.mask(alltrpsF, buffer = 0)

test <- read.mask(data = data.frame(alltrpsF))

mask.buff.F <- Create.buffmask(sigma = 2000, mask = mask)
mask.buff.M <- Create.buffmask(sigma = 5000, mask = mask)

n.reps <- 5
grid.traps.F <- Proposed.traps(alltraps = alltrpsF, nT = 80, D = NULL, sigma = sigmaF, lambda0F,  n.reps = 10, grid = TRUE)

plot(mask.buff.F[[1]])
plot(grid.traps.F[[1]],add = T)

scenF.gridsF <- make.scenarios (trapsindex = 1:n.reps, detectfn = 'HHN', D = DF, lambda0 = lambda0F, 
                                sigma = sigmaF, noccasions = 1)
scenM.gridsF <- make.scenarios (trapsindex = 1:n.reps, detectfn = 'HHN', D = DM, lambda0 = lambda0M, 
                                sigma = sigmaM, noccasions = 1)

nrm.scenF.gridsF <- scenarioSummary(scenF.gridsF, traps = grid.traps.F, mask = mask.buff.F)
nrm.scenM.gridsF <- scenarioSummary(scenM.gridsF, traps = grid.traps.F, mask = mask.buff.M)

###############M
alltrpsM <- make.grid(nx = 50, ny = 50, origin = c(1000,1000), spacing = 2/3 *sigmaM, detector = "count")

grid.traps.M <- Proposed.traps(mask = NULL, alltraps = alltrpsM, nT = 80, D = NULL, sigma = sigmaM, lambda0M,  n.reps = n.reps, grid = TRUE)

scenF.gridsM <- make.scenarios (trapsindex = 1:n.reps, detectfn = 'HHN', D = DF, lambda0 = lambda0F, 
                                sigma = sigmaF, noccasions = 1)
scenM.gridsM <- make.scenarios (trapsindex = 1:n.reps, detectfn = 'HHN', D = DM, lambda0 = lambda0M, 
                                sigma = sigmaM, noccasions = 1)

nrm.scenF.gridsF <- scenarioSummary(scenF.gridsM, traps = grid.traps.M, mask = mask.buff.F)
nrm.scenM.gridsF <- scenarioSummary(scenM.gridsM, traps = grid.traps.M, mask = mask.buff.M)

test <- Proposed.traps(mask = mask.buff.M, alltraps = alltrpsM, nT = 80, D = NULL, sigma = sigmaM, lambda0M,  n.reps = 5, grid = TRUE)

