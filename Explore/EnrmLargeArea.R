#Nov 2023, exploring Enrm summaries for grid and GA designs in a large regular area
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
MDFactor = 1.25  #male to female D

sigmaF = 2500 ; sigmaM = MsigmaFactor * sigmaF ; sigma <- c(sigmaF, sigmaM)
DF = 0.0007 ; DM = MDFactor * DF ; D <- c(DF, DM)
lambda0F = 4 ; lambda0M = 3 ; lambda0 <- c(lambda0F,lambda0M) 

nT <- 40

#create large regular area using approx 30 *M  sigma
#core unbuffered area 30 sigma
#mustn't be a mask object so that the sf functions work
#setting it up to have 0,0 in the centre so in sync with traps for GA

area = st_polygon(  list(
  cbind(c(-75000,75000,75000,-75000,-75000), c(-75000,-75000,75000,75000,-75000 ))
)
)

my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"
test <- st_sf(area, crs = my_crs) 

#use make.grid and make.mask to construct raster mask (and trap locations for GA)
#make.grid places traps at origin (default 0) and adds spacing
#make.mask creates a mask with traps at midpoints
#mask coords are for mid points of cells
#here using the male sigma of 5K
malesig <- 5000
traplocs = make.grid(nx = 45, ny = 45, spacing = 2*malesig/3, detector = "count")
msk <- make.mask(traplocs, buffer = 10000, spacing = malesig/2)

plot(msk, axes = T, border = 0)
plot(traplocs, add = T)

test <- secr::raster(msk)
test_df <- data.frame(coordinates(test))
mask <- read.mask(data = test_df)

#spatial functions
area2 <- st_as_polygon(x = test_df, coords = c(1,2), crs = my_crs) 

x1 <- c(attr(mask,'boundingbox')$x,attr(mask,'boundingbox')$x[1])
y1 <- c(attr(mask,'boundingbox')$y, attr(mask,'boundingbox')$y[1])

test <- st_polygon(list(cbind(x1,y1)))

my_grid <- st_make_grid(area2, cellsize = 100, what = "centers") %>%
  st_intersection(area2)
grid_traps <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)

######################################################################
#single sex designs using supplied grid spacing
#for grids just need to pass the polygon and not the traps obj
#here just using a simple large polygon
grid.F <- Proposed.traps(poly = area, alltraps = NULL, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                         lambda0 = c(lambda0F, lambda0M), buff.sigma = NULL, grid.spacing = sigmaF, 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.M <- Proposed.traps(poly = area, alltraps = NULL, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                         lambda0 = c(lambda0F, lambda0M), buff.sigma = NULL, grid.spacing = sigmaM, 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

grid.Avg <- Proposed.traps(poly = area, alltraps = NULL, D = c(DF, DM), sigma = mean(c(sigmaF,sigmaM)), 
                         lambda0 = c(lambda0F, lambda0M), buff.sigma = NULL, grid.spacing = mean(c(sigmaF,sigmaM)), 
                         criterion = 4, n.reps = 1, grid = TRUE, nT = nT)
plot(area)
plot(grid.F[[1]], add = T)
plot(grid.M[[1]], add = T)
plot(grid.Avg[[1]], add = T)

#GA designs still use a polygon for the regular area
#fn assumes a non-regular area if grid = FALSE and uses type = "polygon" with make.mask
#for the large regular area I will pass an actual buffered mask
#if 2 sets of det pars passed it will use 2 pop approach
#max = F uses the weights for the two Enrm's, T will focus on the minimum of the two

#issues:
#currently an error when using make.mask with a sf obj (only has coords)

#potential locs needed for GA designs
trap.locsF <- make.grid(5,5, spacing = 2/3 * sigmaF, detector = "count")
mask.buffF <- make.mask(trap.locsF, buffer = 25000)

trap.locsM <- make.grid(5,5, spacing = 2/3 * sigmaM, detector = "count")
mask.buffM <- make.mask(trap.locsM, buffer = 25000)

GA.F <- Proposed.traps(poly = area, alltraps = trap.locsF, mask.buff = NULL, D = DF, sigma = sigmaF, 
                         lambda0 = lambda0F, sigma.buff = sigmaF, grid.spacing = NULL, 
                         criterion = 5, n.reps = 1, grid = FALSE, nT = nT)

GA.M <- Proposed.traps(poly = area, alltraps = trap.locsM, mask.buff = mask.buffM, D = DM, sigma = sigmaM, 
                          lambda0 = lambda0M, sigma.buff = sigmaM, grid.spacing = NULL, 
                          criterion = 4, n.reps = 1, grid = FALSE, nT = nT)

#for 2 sex design need to choose what resolution for trap locations
GA.Both.F <- Proposed.traps(poly = area, alltraps = trap.locsF, mask.buff = mask.buffF, D = c(DF, DM), sigma = c(sigmaF,sigmaM), 
                       lambda0 = c(lambda0F, lambda0M), sigma.buff = sigmaF, grid.spacing = NULL, 
                       criterion = 4, n.reps = 1, grid = FALSE, nT = nT)

plot(area)
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





