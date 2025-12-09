#Aug 2025, simulations to evaluate grid design for saigma of 200m with 1600 spacing

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
    sims.group <- fit.models(raw.data, extractfn = predict, fit.args = list(detectfn = 'HHN'),
                             maskset = masklist, fit = TRUE, ncores = numcores, byscenario = F)
  } else {
    sims <- run.scenarios(nreps, scenarios = scenario.df, trapset = traplist, fit = TRUE, extractfn = predict,
                                fit.args = list(detectfn = 'HHN'),
                                maskset = masklist, ncores = numcores, byscenario = F)
  }
  return(list("Data" = raw.data, "Sim results" = sims))
#  return(list("Data summary stats" = data.stats, "CHs" = raw.data, "Sim results" = sims.group))
}

###################################################################
load("SCRObjs.RData")

load("GridDesigns.RData") 

mask <- res.objs[[1]]

#Set up scenario
sigma = 200 
L0 <- 2
D <- 0.05
nT <- 40 ; nreps <- 50
  
scen <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = D, lambda0 = L0, 
                        sigma = sigma, noccasions = 1)

Grid.1600.test <- Sim.designs(scenario.df = scen, nreps = nreps, traplist = grid.designs$`Avg sigma`[[1]], masklist = mask, raw = TRUE, numcores = 20)

save.image("Grid1600Test.RData")

###################################################################
