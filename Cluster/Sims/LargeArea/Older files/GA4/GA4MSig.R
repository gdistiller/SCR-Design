#Sept 2024, simulations to evaluate different proposed designs
#this job is for a GA4 design with M sigma
#rerunning 26 Sept after discovering wrong DF

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
Sim.designs <- function(scenario.df, nreps = 1, traplist, masklist, raw = TRUE, numcores){
  data.stats <- NULL ; raw.data <- NULL
  if (raw==TRUE){
    data.stats <- run.scenarios(nrepl = nreps, scenarios = scenario.df, trapset = traplist, maskset = masklist, extractfn = NULL, fit = FALSE)
    raw.data <- run.scenarios(nrepl = nreps, scenarios = scenario.df, trapset = traplist, maskset = masklist, extractfn = identity, fit = FALSE)
    sims.group <- fit.models(raw.data, extractfn = predict, fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group"),
                             maskset = masklist, fit = TRUE, ncores = numcores)
  } else {
    sims.group <- run.scenarios(nreps, scenarios = scenario.df, trapset = traplist, fit = TRUE, extractfn = predict,
                                fit.args = list(detectfn = 'HHN', model = list(D~g, sigma~g, lambda0~g), groups = "group"),
                                maskset = masklist, ncores = numcores)
  }
  return(list("Data summary stats" = data.stats, "CHs" = raw.data, "Sim results" = sims.group))
}

###################################################################
load("LargeSCRObjs.RData")
load("GADesigns.RData")

mask <- scr.objs[[1]]

#Use CLT type values with 40 traps
MsigmaFactor = 3 #male to female sigma
MDFactor = 0.875  #male to female D

sigmaF = 2000 ; sigmaM = MsigmaFactor * sigmaF ; sigma <- c(sigmaF, sigmaM)
DF = 0.0008 ; DM = MDFactor * DF ; D <- c(DF, DM)
lambda0F = 6 ; lambda0M = 3 ; lambda0 <- c(lambda0F,lambda0M) 
nT <- 40 ; nreps = 500

scen <- make.scenarios (trapsindex = 1:nreps, detectfn = 'HHN', D = DF, lambda0 = lambda0F, 
                        sigma = sigmaF, noccasions = 1, groups = c('F','M'))

male <- scen$group == 'M'
scen$D[male] <- DM
scen$lambda0[male] <- lambda0M
scen$sigma[male] <- sigmaM

GA4.MSigma.results <- Sim.designs(scenario.df = scen, nreps = 1, traplist = GA.designs$`M4 Sigma`, masklist = mask, raw = TRUE, numcores = 10)

save.image("GA4MSigma.RData")

###################################################################







