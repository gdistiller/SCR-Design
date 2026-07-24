#July 2026, simulations to evaluate 2 stage designs for 15 fold scen, for 120 traps
#rerunning as versionb c in case there is an issue with the current jobs
#there is a (from ngen = 300) and b but not sure if the same R script will cause problems

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

#sim data
#uses secrdesign to simulate data
Sim.data <- function(scenario.df, nreps = 1, traplist, masklist){
  raw.data <- run.scenarios(nrepl = nreps, scenarios = scenario.df, trapset = traplist, maskset = masklist, extractfn = identity, fit = FALSE, seed = 2025)
  return("Data" = raw.data)
}

#fit models without passing start values
#Note: change in construction of CH, due to only one trap being used for LW and how secrdesign formats the data obj
Fit.models <- function(Data.obj, nrep, numcores = 4){
  ests <- matrix(NA, nrow = nrep, ncol = 12, 
                 dimnames = list(NULL, c("D1", "D1.se", "L01", "L01.se", "Sig1", "Sig1.se", "D2", "D2.se", "L02", "L02.se", "Sig2", "Sig2.se")))
  for (i in 1:nrep){
    CH <- Data.obj$output[[i]][[1]]
    mod <- secr.fit(CH, detectfn = 'HHN', mask = mask, model = list(D~g, sigma~g, lambda0~g), groups = "group", ncores = numcores, trace = FALSE)
    temp1 <- unlist(predict(mod)[[1]][,2:3])
    temp2 <- unlist(predict(mod)[[2]][,2:3])
    ests[i,] <- c(temp1[c(1,4,2,5,3,6)], temp2[c(1,4,2,5,3,6)])
  }
  return(ests)
}

###################################################################

load("SCRObjs.RData")
load("GA2StageDesignsb.RData") #XXb is for designs from ngen = 500

mask <- SCR.objs$Full[[1]]

#Set values for both strata
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor ; L0 <- c(L01,L02)
D1 <- 0.05 ; D2 <- D1 / DiffFactor ; D <- c(D1, D2)
nT <- 120 ; nreps <- 500 ; cores <- 40

#only 1 trap index, multiple rows with 1 trap and maskindex throws an error
scen <- make.scenarios (trapsindex = 1:nreps, detectfn = 'HHN', D = D1, lambda0 = L01, 
                        sigma = sigma1, noccasions = 1, groups = c('S1','S2'))

S2 <- scen$group == 'S2'
scen$D[S2] <- D2
scen$lambda0[S2] <- L02
scen$sigma[S2] <- sigma2

TwoStage.120.Data <- Sim.data(scenario.df = scen, nreps = 1, traplist = GA2StageDesigns$`120 traps`, masklist = mask)
TwoStage.120.results <- Fit.models(Data.obj = TwoStage.120.Data, nrep = nreps, numcores = cores)
save.image("TwoStageResults120c.RData")


###################################################################

