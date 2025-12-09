#July 2025, simulations to evaluate grid designs for 15 fold scen
#pass secrdesign a list of traps so can't use foreach
#using 20 cores
#grids use spacing of either 800m (4 * small sigma or about 1/4 large sigma) 
#or 1600 m (avg sigma)
#modified in August to add start = NULL to remove effect of starting at true values (versionX2.RData)

rm(list=ls())

#################### INITIALIZE MPI ENVIRONMENT ####################
library(secr)
library(secrdesign)

# In case R exits unexpectedly, have it automatically clean up resources taken up by Rmpi (slaves, memory, etc...) 
# This just provides some guarantee that if R is exited before closing the slaves, and without cleaning up the MPI environment, it will do this for you.
.Last = function(){
  if (is.loaded("mpi_initialize"))
  {
    if (mpi.comm.size(1) > 0)
    {
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}

#sim function
#note that one can also get the summary stats from the CHs so not necessary to save both 
Sim.designs <- function(scenario.df, nreps = 1, traplist, masklist, raw = TRUE, numcores){
  data.stats <- NULL ; raw.data <- NULL
  if (raw==TRUE){
#    data.stats <- run.scenarios(nrepl = nreps, scenarios = scenario.df, trapset = traplist, maskset = masklist, extractfn = NULL, fit = FALSE)
    raw.data <- run.scenarios(nrepl = nreps, scenarios = scenario.df, trapset = traplist, maskset = masklist, extractfn = identity, fit = FALSE)
    sims.group <- fit.models(raw.data, extractfn = predict, fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group", start = NULL),
                             maskset = masklist, fit = TRUE, ncores = numcores, byscenario = F)
  } else {
    sims.group <- run.scenarios(nreps, scenarios = scenario.df, trapset = traplist, fit = TRUE, extractfn = predict,
                                fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group", start = NULL),
                                maskset = masklist, ncores = numcores, byscenario = F)
  }
  return(list("Data" = raw.data, "Sim results" = sims.group))
#  return(list("Data summary stats" = data.stats, "CHs" = raw.data, "Sim results" = sims.group))
}

###################################################################
load("SCRObjs.RData")

load("GridDesigns.RData") 

mask <- res.objs[[1]]

#Set values for both strata
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor ; L0 <- c(L01,L02)
D1 <- 0.05 ; D2 <- D1 / DiffFactor ; D <- c(D1, D2)
nT <- 40 ; nreps <- 100
  
scen <- make.scenarios (trapsindex = 1:nreps, detectfn = 'HHN', D = D1, lambda0 = L01, 
                        sigma = sigma1, noccasions = 1, groups = c('S1','S2'))

S2 <- scen$group == 'S2'
scen$D[S2] <- D2
scen$lambda0[S2] <- L02
scen$sigma[S2] <- sigma2

Grid.800.results <- Sim.designs(scenario.df = scen, nreps = 1, traplist = grid.designs$`800 m`, masklist = mask, raw = TRUE, numcores = 20)
save.image("GridResults240.RData")

Grid.1600.results <- Sim.designs(scenario.df = scen, nreps = 1, traplist = grid.designs$`Avg sigma`, masklist = mask, raw = TRUE, numcores = 20)
save.image("GridResults240.RData")

###################################################################
