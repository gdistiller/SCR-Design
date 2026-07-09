#Aug 2024, refitting models to CLT
#fitting SCR models to CLT data
#includes corrected traps and the few NA sex individuals
##this version uses the complex mask without a buffer and 1 occ
rm(list=ls())

#################### INITIALIZE MPI ENVIRONMENT ####################
library(secr)

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

###################################################################
#load data, traps and mask

setwd("/home/gdistiller/Design/CLT/")

load("CLTAreas.RData")
load("CLTCapthists.RData")

#both sessions together
fit.all.null <- secr.fit(CH.red, mask = CLTcomp.nobuff, hcov = "Sex",list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.all.sex <- secr.fit(CH.red, mask = CLTcomp.nobuff, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

save.image("CLTResultsRed.RData")

#seperately using MS obj
#note need to specify pmix ~ session to est pmix per area (b)
fit.MS.NS.null <- secr.fit(MS.red, mask = CLTcomp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session, sigma ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexa <- secr.fit(MS.red, mask = CLTcomp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session + h2, sigma ~ session + h2), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexb <- secr.fit(MS.red, mask = CLTcomp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session + h2, sigma ~ session + h2, pmix ~ session + h2), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexc <- secr.fit(MS.red, mask = CLTcomp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2, pmix ~ session + h2), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexd <- secr.fit(MS.red, mask = CLTcomp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2, pmix ~ session*h2), detectfn = 'HHN', trace = FALSE, ncores = 1)

save.image("CLTResultsRed2.RData")

###################################################################
