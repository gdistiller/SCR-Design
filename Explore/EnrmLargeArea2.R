#August 2024
#create large regular area using approx 30 * M  sigma
#core unbuffered area 30 sigma
#creating area as both polygon and raster object
#polygon required for randomised grid designs, raster -> traps / mask required for GA designs

library(ggplot2)
library(secr)
library(secrdesign)
library(kofnGA)
library(dplyr)
library(tidyr)

###############################
#using sf, used for vector data
###############################

library(sf)

#create sf object using st_sf()
#need to provide two elements, a df with attributes and a simple feature geometry list-col
#1st create simple feature geometries sfg, then use st_sfc(), and then st_sf to combine with df
p1 <- rbind(c(0, 0), c(150000, 0), c(150000, 150000), c(0, 150000), c(0,0))
area.poly <- st_polygon(list(p1)) #creates simple feature geometries sfg
my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"
area.sfc <- st_sfc(area.poly, crs = my_crs) #creates simple feature geometry list-column

df <- data.frame(v1 = c("Polygon mask"))  #df needed to create sf object
area <- st_sf(df, geometry = area.sfc) 
head(area)

plot(area, axes = T)

##################################################
#using terra, used for both vector and raster data
##################################################

library(terra)

#create raster for mask and for trap locations
mask.r <- rast(ncol = 120, nrow = 120, xmin = 0, xmax = 150000, ymin = 0, ymax = 150000)  #ncol/row determine the resolution
traplocs.r <- rast(ncol = 60, nrow = 60, xmin = 30000, xmax = 120000, ymin = 30000, ymax = 120000) #res of 1.5K

nrow(mask.r)
dim(mask.r)
ncell(mask.r)

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

plot(mask, axes = T)
plot(trap.locs, add = T)

#try out 40 traps

D4 = 0.0004 ; L0 = 3 ; sig = 3000

nT <- 40

######
#Grid
######

#single sex designs using supplied grid spacing
#for grids just need to pass the polygon (for trap locs)
locs <- st_as_sf(traplocs.v)
msk <- st_as_sf(mask.v)

grid.15 <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                         lambda0 = L0, sigma.buff = NULL, grid.spacing = 1500, 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.20 <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                          lambda0 = L0, sigma.buff = NULL, grid.spacing = 2000, 
                          criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

plot(mask.v)
plot(grid.15[[1]], add = T)

######
#GA
######
#GA designs still use a polygon for the regular area
#fn assumes a non-regular area if grid = FALSE and uses type = "polygon" with make.mask
#for the large regular area I will pass an actual buffered mask
#if 2 sets of det pars passed it will use 2 pop approach
#max = F uses the weights for the two Enrm's, T will focus on the minimum of the two

GA.F4 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D4, sigma = sig, 
                       lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                       criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 20)
plot(mask.v)
plot(GA.F4[[1]], add = T)

GA.F5 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask,  D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 5, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 20)
plot(mask.v)
plot(GA.F5[[1]], add = T)

#Enrm numbers
#scenarios with F / M SCR parameters
scenF <- make.scenarios (trapsindex = 1:3, detectfn = 'HHN', D = 0.0004, lambda0 = 3, 
                         sigma = 3000, noccasions = 1)

traplist <- list(Grid.40 = grid.15[[1]], GA.fill = GA.F4[[1]], GA.cluster = GA.F5[[1]])

scenarioSummary(scenF, traps = traplist, mask = mask)

test.data <- run.scenarios(10, scenF, traplist, mask, fit = FALSE, extract = identity, 
                           det.args = list(savepopn = TRUE), fit.args = list(detectfn = 'HHN'))
summary(test.data)

#can save data and simulated pop
pop1 <- attr(test.data$output[[1]][[1]],"popn")
pop2 <- attr(test.data$output[[2]][[1]],"popn")

plot(mask, dots = FALSE)
plot(traplist[[1]], detpar=list(col="green", pch = 15), add = TRUE)
plot(pop1, frame = FALSE, add = TRUE, col = 'blue', pch = 16, cex = 0.6)

str(test.data$output[[1]][[1]])
test.data$output[[1]][[1]][,1,]
apply(test.data$output[[1]][[1]][,1,],1,sum)

test.fit.40 <- run.scenarios(20, scenF, traplist, mask, fit = TRUE, fit.args = list(detectfn = 'HHN'))

#extract sim results. select.stats can give rel bias, rel se, abs deviation, coverage
find.param(test.fit)

stats1 <- select.stats(test.fit, parameter = "D", statistics = c("estimate","lcl", "ucl", "RB", "RSE", "COV"))
lapply(stats1$output, head, 4)

#can filter rogue values, need to pass selected stats obj

x <- validate (x = stats1, test = "estimate", validrange = c(0, 0.003), targets = "all")

summary (x, dec = 5, fields = c("n", "mean", "se"), alpha = 0.05,
         type = c("list", "dataframe", "array"))

#Applying the ‘rms’ field to the absolute deviation of an estimate (ERR) provides the RMSE.

#plot
par(mfrow = c(3,1))
plot(stats1, type = "hist", statistic = "estimate")
plot(stats1, type = "CI")

parallel::detectCores()

for (i in 1:3){
  boxplot(x$output[[i]][[1]], ylim=c(0,0.0006))
  abline(h = 0.0004, col = "red", lwd = 2)
}

##########
#continue to test functions
#Now GA designs with no mask, buffered mask created from traplocs with type = "trapbuffer"
#for irregular area, if want custom area must provide mask.buff

GA.F4.test <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = NULL, D = D, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 4, n.reps = 1, grid = FALSE, nT = nT)

plot(mask.v)
plot(GA.F4.test[[1]], add = T)

#############
#Now with an irregular area 

load("CLTAreas.RData")

#Need CLTsimp.buff mask and CLTsimp_trap_locs_Avg traps obj, and also the polygon
#This polygon covers the whole space for trap locs

plot(CLTsimp.buff)
plot(CLTsimp_trap_locs_Avg, add = T)

#grid
grid.15 <- Proposed.traps(poly = CLTSimple$geometry, alltraps = NULL, D = D4, sigma = sig, 
                          lambda0 = L0, sigma.buff = NULL, grid.spacing = 1500, 
                          criterion = 4, n.reps = 1, grid = TRUE, nT = nT)
plot(CLTsimp.buff)
plot(grid.15[[1]], add = T)

#GA with mask
GA.F4 <- Proposed.traps(poly = NULL, alltraps = CLTsimp_trap_locs_Avg, mask.buff = CLTsimp.buff, D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 20)
plot(CLTsimp.buff)
plot(GA.F4[[1]], add = T)

GA.F5 <- Proposed.traps(poly = NULL, alltraps = CLTsimp_trap_locs_Avg, mask.buff = CLTsimp.buff, D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 5, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 20)

plot(CLTsimp.buff)
plot(GA.F5[[1]], add = T)

#GA without mask
GA.F4 <- Proposed.traps(poly = NULL, alltraps = CLTsimp_trap_locs_Avg, mask.buff = NULL, D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 20)
plot(CLTsimp.buff)
plot(GA.F4[[1]], add = T)

GA.F5 <- Proposed.traps(poly = NULL, alltraps = CLTsimp_trap_locs_Avg, mask.buff = NULL, D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 5, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 20)

plot(CLTsimp.buff)
plot(GA.F5[[1]], add = T)

############################################
#continue with large area
#look at spacing
############################################
#L grid sizes, 1500 spacing
grid.15L <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                           lambda0 = L0, sigma.buff = NULL, grid.spacing = 1500, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 80)

#L grid sizes, 2000 spacing
grid.20L <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                           lambda0 = L0, sigma.buff = NULL, grid.spacing = 2000, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 80)

#L grid sizes, 3000 spacing
grid.30L <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                           lambda0 = L0, sigma.buff = NULL, grid.spacing = 3000, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 80)

#L grid sizes, 4000 spacing
grid.40L <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                           lambda0 = L0, sigma.buff = NULL, grid.spacing = 4000, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 80)

#L grid sizes, 5000 spacing
grid.50L <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                           lambda0 = L0, sigma.buff = NULL, grid.spacing = 5000, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 80)

#L grid sizes, 5000 spacing
grid.60L <- Proposed.traps(poly = locs, alltraps = NULL, D = D4, sigma = sig, 
                           lambda0 = L0, sigma.buff = NULL, grid.spacing = 6000, 
                           criterion = 4, n.reps = 1, grid = TRUE, nT = 80)

scen.grids.L <- make.scenarios (trapsindex = 1:6, detectfn = 'HHN', D = D4, lambda0 = L0, 
                         sigma = sig, noccasions = 1)

traplist <- list(Grid15.80 = grid.15L[[1]], Grid20.80 = grid.20L[[1]], Grid30.80 = grid.30L[[1]], 
                 Grid40.80 = grid.40L[[1]], Grid50.80 = grid.50L[[1]], Grid60.80 = grid.60L[[1]])

scenarioSummary(scen.grids.L, traps = traplist, mask = mask)

grids.L.data <- run.scenarios(50, scen.grids.L, traplist, mask, fit = FALSE, fit.args = list(detectfn = 'HHN'))
summary(grids.L.data)

grids.L.fit <- run.scenarios(10, scen.grids.L, traplist, mask, fit = TRUE, fit.args = list(detectfn = 'HHN'))

#extract sim results. select.stats can give rel bias, rel se, abs deviation, coverage
grid.stats1 <- select.stats(grids.L.fit, parameter = "D", statistics = c("estimate","lcl", "ucl", "RB", "RSE", "COV"))
x <- validate (x = grid.stats1, test = "estimate", validrange = c(0, 0.003), targets = "all")

summary (x, dec = 5, fields = c("n", "mean", "se"), alpha = 0.05,
         type = c("list", "dataframe", "array"))

#plot
par(mfrow = c(2,2))
plot(grid.stats1, type = "hist", statistic = "estimate")
plot(grid.stats1, type = "CI")

for (i in 1:3){
  boxplot(x$output[[i]][[1]], ylim=c(0,0.0006))
  abline(h = 0.0004, col = "red", lwd = 2)
}

#######################
#explore GA designs, data summaries

GA.F4.traplist <- vector("list", 20)

for (i in 1:20){
  GA.F4 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D4, sigma = sig, 
                          lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                          criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 50)
  GA.F4.traplist[[i]] <- GA.F4[[1]] 
}


GA.F5.traplist <- vector("list", 20)

for (i in 1:20){
  GA.F5 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D4, sigma = sig, 
                          lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                          criterion = 5, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 50)
  GA.F5.traplist[[i]] <- GA.F5[[1]] 
}

GA.traps <- list(Crit4 = GA.F4.traplist, Crit5 = GA.F5.traplist)

save(GA.traps, file = "GAtraps.RData")

scen.GA.F4.40 <- make.scenarios (trapsindex = 1:20, detectfn = 'HHN', D = D4, lambda0 = L0, 
                                sigma = sig, noccasions = 1)

scenarioSummary(scen.GA.F4.40, traps = GA.F4.traplist, mask = mask)

GA.F4.data <- run.scenarios(50, scen.GA.F4.40, GA.F4.traplist, mask, fit = FALSE, fit.args = list(detectfn = 'HHN'))
summary(GA.F4.data)

scen.GA.F5.40 <- make.scenarios (trapsindex = 1:20, detectfn = 'HHN', D = D4, lambda0 = L0, 
                                 sigma = sig, noccasions = 1)

scenarioSummary(scen.GA.F5.40, traps = GA.F5.traplist, mask = mask)

GA.F5.data <- run.scenarios(50, scen.GA.F5.40, GA.F5.traplist, mask, fit = FALSE, fit.args = list(detectfn = 'HHN'))
summary(GA.F5.data)

#fit
GA.F4.fit <- run.scenarios(1, scen.GA.F4.40, GA.F4.traplist, mask, fit = TRUE, fit.args = list(detectfn = 'HHN'), ncores = 10)
summary(GA.F4.fit)

Dests <- vector("numeric", 20)
for (i in 1:20){
  Dests[[i]] <- GA.F4.fit$output[[i]][[1]][1,2]
}

boxplot(Dests)
summary(Dests)

(mean(Dests)-0.0004) / 0.0004

GA.F5.fit <- run.scenarios(1, scen.GA.F5.40, GA.F5.traplist, mask, fit = TRUE, fit.args = list(detectfn = 'HHN'), ncores = 10)
summary(GA.F5.fit)

Dests2 <- vector("numeric", 20)
for (i in 1:20){
  Dests2[[i]] <- GA.F5.fit$output[[i]][[1]][1,2]
}

boxplot(Dests2)
abline(h = 0.0004)
summary(Dests2)

(mean(Dests2)-0.0004) / 0.0004

#rerunning this time with one realisation from the GA design
#seems like less variance as one would expect 
scen.GA.40 <- make.scenarios (trapsindex = 1:2, detectfn = 'HHN', D = D4, lambda0 = L0, 
                                 sigma = sig, noccasions = 1)

scenarioSummary(scen.GA.40, traps = list(Crit4 = GA.F4.traplist[[1]], Crit5 = GA.F5.traplist[[1]]), mask = mask)
GA.data <- run.scenarios(20, scen.GA.40, list(Crit4 = GA.F4.traplist[[1]], Crit5 = GA.F5.traplist[[1]]), mask, fit = FALSE, fit.args = list(detectfn = 'HHN'))
GA.fit <- run.scenarios(20, scen.GA.40, list(Crit4 = GA.F4.traplist[[1]], Crit5 = GA.F5.traplist[[1]]), mask, fit = TRUE, fit.args = list(detectfn = 'HHN'), ncores = 10)


#################################################################################
###two groups
#set SCR parameters
MsigmaFactor = 2 #male to female sigma
MDFactor = 1.25  #male to female D

sigmaF = 2500 ; sigmaM = MsigmaFactor * sigmaF ; sigma <- c(sigmaF, sigmaM)
DF = 0.0007 ; DM = MDFactor * DF ; D <- c(DF, DM)
lambda0F = 4 ; lambda0M = 3 ; lambda0 <- c(lambda0F,lambda0M) 






#for 2 sex design need to choose what resolution for trap locations
GA.Both.F <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                            lambda0 = c(lambda0F, lambda0M), sigma.buff = sigmaF, grid.spacing = NULL, 
                            criterion = 4, n.reps = 1, grid = FALSE, nT = nT)


#try when a buffered mask is not passed
GA.F <- Proposed.traps(poly = mask.v, alltraps = trap.locs, mask.buff = NULL, D = DF, sigma = sigmaF, 
                       lambda0 = lambda0F, sigma.buff = sigmaF, grid.spacing = NULL, 
                       criterion = 4, n.reps = 1, grid = FALSE, nT = nT)



plot(area)
plot(GA.F[[1]], add = T)
plot(GA.M[[1]], add = T)


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



#estimateSummary useful for groups
estimateSummary(x, 'D',  format = 'data.frame', cols = c(1,2,6:8))
estimateSummary(test.fit, 'lambda0',  format = 'data.frame',cols = c(1,2,6:8))
estimateSummary(test.fit, 'sigma',  format = 'data.frame',cols = c(1,2,6:8))



########################################################################
#non-U D
#need mask with density covariate
#pop.args must include model2D = "IHP" and D = "XX" where XX is the covariate
#can visualise simulate dpopns with savepopn = TRUE in det.args
#use popindex to compare different pop.args (can have different Ds as different mask covariates)


