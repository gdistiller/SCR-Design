#Proposed designs for CLT area
#redone below with revised mask and estimates

load("LargeSCRObjs.RData")

#set up scenario
MsigmaFactor = 3 #male to female sigma
MDFactor = 0.875  #male to female D

sigmaF = 2000 ; sigmaM = MsigmaFactor * sigmaF ; sigma <- c(sigmaF, sigmaM)
DF = 0.00008 ; DM = MDFactor * DF ; D <- c(DF, DM)
lambda0F = 6 ; lambda0M = 3 ; lambda0 <- c(lambda0F,lambda0M) 
nT <- 40

Scen.F <- make.scenarios (trapsindex = 1, D = DF, lambda0 = lambda0F, sigma = sigmaF, noccasions = 1, detectfn = 14) 
Scen.M <- make.scenarios (trapsindex = 1, D = DM, lambda0 = lambda0M, sigma = sigmaM, noccasions = 1, detectfn = 14) 

#load polygon, trap locs and mask
load("CLT/SepAreas.RData")

CLTgrid.2sigF <- Proposed.traps(poly = Simp_utm, alltraps = NULL, D = NULL, sigma = NULL, 
                             lambda0 = NULL, sigma.buff = NULL, grid.spacing = sigmaF, 
                             criterion = 4, n.reps = 1, grid = TRUE, nT = nT)

plot(st_geometry(Simp_utm))
plot(CLTgrid.2sigF[[1]], add = T)

#for GA design
GA4test <- Proposed.traps(poly = NULL, alltraps = CLTsimp_trap_locs, mask.buff = CLTsimp.nobuff, D = DF, sigma = sigmaF, 
               lambda0 = lambda0F, sigma.buff = sigmaF, grid.spacing = NULL, 
               criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 50)

GA5test <- Proposed.traps(poly = NULL, alltraps = CLTsimp_trap_locs, mask.buff = CLTsimp.nobuff, D = DF, sigma = sigmaF, 
                          lambda0 = lambda0F, sigma.buff = sigmaF, grid.spacing = NULL, 
                          criterion = 5, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 50)


##############################################################################
#creating grid designs for CLT areas
load("CLT/SepAreas.RData")
nreps = 500

#South
CLT.S.Grids.F3K.list <- vector("list", nreps)
CLT.S.Grids.M5K.list <- vector("list", nreps)

for (i in 1:nreps){
  
CLT.S.Grids.F3K.list[[i]] <- Proposed.traps(poly = CLTmsks$`S Polygon SF obj`, alltraps = CLTmsks$`Trap locs for South`, mask.buff = CLTmsks$`Mask of South`, D = 0.00009, sigma = 2000, 
                 lambda0 = 6, sigma.buff = NULL, grid.spacing = 3000, criterion = 4, 
                 n.reps = 1, grid = T, nT = 40, GA.ngen = 50)

CLT.S.Grids.M5K.list[[i]] <- Proposed.traps(poly = CLTmsks$`S Polygon SF obj`, alltraps = CLTmsks$`Trap locs for South`, mask.buff = CLTmsks$`Mask of South`, D = 0.00006, sigma = 6500, 
                                 lambda0 = 3.25, sigma.buff = NULL, grid.spacing = 5000, criterion = 4, 
                                 n.reps = 1, grid = T, nT = 40, GA.ngen = 50)

}

CLT.S.grids <- list("CLT S Grid 3K" = CLT.S.Grids.F3K.list, "CLT S Grid 5K" = CLT.S.Grids.M5K.list) 

#North
CLT.N.Grids.F3K.list <- vector("list", nreps)
CLT.N.Grids.M5K.list <- vector("list", nreps)

for (i in 1:nreps){
  
  CLT.N.Grids.F3K.list[[i]] <- Proposed.traps(poly = CLTmsks$`N Polygon SF obj`, alltraps = CLTmsks$`Trap locs for North`, mask.buff = CLTmsks$`Mask of North`, D = 0.00008, sigma = 2500, 
                                              lambda0 = 0.03, sigma.buff = NULL, grid.spacing = 3000, criterion = 4, 
                                              n.reps = 1, grid = T, nT = 40, GA.ngen = 50)
  
  CLT.N.Grids.M5K.list[[i]] <- Proposed.traps(poly = CLTmsks$`N Polygon SF obj`, alltraps = CLTmsks$`Trap locs for North`, mask.buff = CLTmsks$`Mask of North`, D = 0.0001, sigma = 3700, 
                                              lambda0 = 0.0325, sigma.buff = NULL, grid.spacing = 5000, criterion = 4, 
                                              n.reps = 1, grid = T, nT = 40, GA.ngen = 50)
  
}

CLT.N.grids <- list("CLT N Grid 3K" = CLT.N.Grids.F3K.list, "CLT N Grid 5K" = CLT.N.Grids.M5K.list) 

save(CLT.S.grids,CLT.N.grids, file = "CLTSepGrids.RData")

#GA designs on cluster

