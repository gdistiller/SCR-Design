#May 2026, simulations to evaluate different proposed GA4 designs for 15 fold scen, for 120 traps
#version d is for both strata pars
#designs after increasing ngen
#using my own function to fit models to data simulated by secrdesign, 
#tried to use foreach but can't get nice collation of results
#reverting to rather using numcores
#this version does 500 reps and uses 25 cores

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
  raw.data <- run.scenarios(nrepl = nreps, scenarios = scenario.df, trapset = traplist, maskset = masklist, extractfn = identity, fit = FALSE, seed = 2026)
  return("Data" = raw.data)
}

#fit models without passing start values
Fit.models <- function(Data.obj, nrep, numcores = 4){
  ests <- matrix(NA, nrow = nrep, ncol = 12, 
                 dimnames = list(NULL, c("D1", "D1.se", "L01", "L01.se", "Sig1", "Sig1.se", "D2", "D2.se", "L02", "L02.se", "Sig2", "Sig2.se")))
  for (i in 1:nrep){
    CH <- Data.obj$output[[i]][[1]]
    mod <- secr.fit(CH, detectfn = 'HHN', mask = mask, model = list(D~g, sigma~g, lambda0~g), groups = "group", ncores = numcores, trace = F)
    temp1 <- unlist(predict(mod)[[1]][,2:3])
    temp2 <- unlist(predict(mod)[[2]][,2:3])
    ests[i,] <- c(temp1[c(1,4,2,5,3,6)], temp2[c(1,4,2,5,3,6)])
    print(i)
  }
  return(ests)
}

###################################################################
load("SCRObjs.RData")

load("GADesigns120.RData") 

mask <- SCR.objs$Full[[1]]

#Set values for both strata
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor ; L0 <- c(L01,L02)
D1 <- 0.05 ; D2 <- D1 / DiffFactor ; D <- c(D1, D2)
nT <- 120 ; nreps <- 500 ; cores <- 25
  
scen <- make.scenarios (trapsindex = 1:nreps, detectfn = 'HHN', D = D1, lambda0 = L01, 
                        sigma = sigma1, noccasions = 1, groups = c('S1','S2'))

S2 <- scen$group == 'S2'
scen$D[S2] <- D2
scen$lambda0[S2] <- L02
scen$sigma[S2] <- sigma2

G5.Both.Datad <- Sim.data(scenario.df = scen, nreps = 1, traplist = GA.designs.120$G5Both, masklist = mask)
G5.Both.resultsd <- Fit.models(Data.obj = G5.Both.Datad, nrep = nreps, numcores = cores)
save.image("GA5Results120d.RData")

###################################################################
