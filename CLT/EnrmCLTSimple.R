#Sept 2023, exploring Enrm summaries for grid and GA designs in the CLT simple area
library(secr)
library(secrdesign)
library(kofnGA)
library(raster)
library(dplyr)
library(tidyr)
library(sf)

#################################################################
#use Proposed traps function to get a grid design and a GA design
#################################################################

#set SCR parameters
MsigmaFactor = 2 #male to female sigma
MDFactor = 0.2  #male to female D

sigmaF = 2400 ; sigmaM = MsigmaFactor * sigmaF ; sigma <- c(sigmaF, sigmaM)
DF = 0.0007 ; DM = MDFactor * DF ; D <- c(DF, DM)
lambda0F = 8 ; lambda0M = 6 ; lambda0 <- c(lambda0F,lambda0M) 

nT <- 40

#get CLT simple study area
#includes traplocs and buffered masks for different sigmas
load("CLTAreas.RData")

######################################################################
#single sex designs using 1 sigma grid spacing
#for grids need to pass the polygon and not the traps obj
grid.F <- Proposed.traps(poly = CLTSimple, alltraps = NULL, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                         lambda0 = c(lambda0F, lambda0M), buff.sigma = NULL, grid.spacing = sigmaF, 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.M <- Proposed.traps(poly = CLTSimple, alltraps = NULL, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                         lambda0 = c(lambda0F, lambda0M), buff.sigma = NULL, grid.spacing = sigmaM, 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.Avg <- Proposed.traps(poly = CLTSimple, alltraps = NULL, D = c(DF, DM), sigma = mean(c(sigmaF,sigmaM)), 
                         lambda0 = c(lambda0F, lambda0M), buff.sigma = NULL, grid.spacing = mean(c(sigmaF,sigmaM)), 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

#GA designs with CLT area as shape file passed
#fn assumes a non-regular area if grid = FALSE and uses type = "polygon" with make.mask
#if 2 sets of det pars passed it will use 2 pop approach
#max = F uses the weights for the two Enrm's, T will focus on the minimum of the two
GA.F <- Proposed.traps(poly = CLTSimple, alltraps = CLTsimp_trap_locs_SigmaF, D = DF, sigma = sigmaF, 
                         lambda0 = lambda0F, buff.sigma = sigmaF, grid.spacing = NULL, 
                         criterion = 4, n.reps = 1, grid = FALSE, nT = nT)

GA.M <- Proposed.traps(poly = CLTSimple, alltraps = CLTsimp_trap_locs_SigmaM, D = DM, sigma = sigmaM, 
                          lambda0 = lambda0M, buff.sigma = sigmaM, grid.spacing = NULL, 
                          criterion = 4, n.reps = 1, grid = FALSE, nT = nT)

GA.Both <- Proposed.traps(poly = CLTSimple, alltraps = CLTsimp_trap_locs_SigmaM, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                       lambda0 = c(lambda0F, lambda0M), buff.sigma = sigmaM, grid.spacing = NULL, 
                       criterion = 4, n.reps = 1, grid = FALSE, nT = nT)


plot(st_geometry(CLTSimple), axes = T)
plot(grid.F[[1]], add = T)
plot(grid.M[[1]], add = T)
plot(grid.Avg[[1]], add = T)

plot(GA.F[[1]], add = T)
plot(GA.M[[1]], add = T)

#plot for SEEC research meeting
pdf("ProposedGrid1.pdf", width = 6, height = 9)
plot(st_geometry(CLTSimple), axes = T)
plot(grid.F[[1]], add = T)
dev.off()

pdf("ProposedGA1.pdf", width = 6, height = 9)
plot(st_geometry(CLTSimple), axes = T)
plot(GA.F[[1]], add = T)
dev.off()

pdf("ProposedGABoth.pdf", width = 6, height = 9)
plot(st_geometry(CLTSimple), axes = T)
plot(GA.F[[1]], add = T)
dev.off()

#now get Enrm summaries, one fine integration mask used

#scenarios with F / M SCR parameters
scenF <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = DF, lambda0 = lambda0F, 
                                sigma = sigmaF, noccasions = 1)

scenM <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = DM, lambda0 = lambda0M, 
                                sigma = sigmaM, noccasions = 1)

#results for grid designs
nrm.scenF.gridsF <- scenarioSummary(scenF, traps = grid.F, mask = CLTsimp.nobuff)
nrm.scenM.gridsF <- scenarioSummary(scenM, traps = grid.F, mask = CLTsimp.nobuff)

nrm.scenF.gridsM <- scenarioSummary(scenF, traps = grid.M, mask = CLTsimp.nobuff)
nrm.scenM.gridsM <- scenarioSummary(scenM, traps = grid.M, mask = CLTsimp.nobuff)

#results for GA designs
nrm.scenF.GA.F <- scenarioSummary(scenF, traps = GA.F, mask = CLTsimp.nobuff)
nrm.scenM.GA.F <- scenarioSummary(scenM, traps = GA.F, mask = CLTsimp.nobuff)

nrm.scenF.GA.M <- scenarioSummary(scenF, traps = GA.M, mask = CLTsimp.buff)
nrm.scenM.GA.M <- scenarioSummary(scenM, traps = GA.M, mask = CLTsimp.buff)





