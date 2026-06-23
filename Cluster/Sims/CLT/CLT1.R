#Sept 2023
#fitting SCR models to CLT data
#this with the 3 traps not found renamed to be R5
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

#first fitting to all data as 1 session
fit.null <- secr.fit(CH, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex <- secr.fit(CH, mask = CLTsimp.nobuff, list(D ~ g, lambda0 ~ g, sigma ~ g), groups = "Sex", detectfn = 'HHN')

save.image("CLTResults.RData")

#session 1 fits
fit.null.session1 <- secr.fit(CH.session1, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex.session1 <- secr.fit(CH.session1, mask = CLTsimp.nobuff, list(D ~ g, lambda0 ~ g, sigma ~ g), groups = "Sex", detectfn = 'HHN')

save.image("CLTResults.RData")

#session 2 fits
fit.null.session2 <- secr.fit(CH.session2, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex.session2 <- secr.fit(CH.session2, mask = CLTsimp.nobuff, list(D ~ g, lambda0 ~ g, sigma ~ g), groups = "Sex", detectfn = 'HHN')

save.image("CLTResults.RData")
###################################################################
