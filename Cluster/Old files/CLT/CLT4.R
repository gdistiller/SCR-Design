#Sept 2024
#fitting SCR models to CLT data
#this with the 3 traps not found renamed to be R5
#redone with revised mask from Anita
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

setwd("/home/01376436/Design/CLT/")

load("SepAreas.RData")
load("CLTCapthists.RData")

#first using multiple occasions
#fitting to all data as 1 session
fit.null.both <- secr.fit(CH, mask = CLT_Both_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex.both <- secr.fit(CH, mask = CLT_Both_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

save.image("CLTResults2.RData")

#N fits
fit.null.N <- secr.fit(CH.Northern, mask = CLT_N_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex.N <- secr.fit(CH.Northern, mask = CLT_N_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2),  detectfn = 'HHN')

save.image("CLTResults2.RData")

#S fits
fit.null.S <- secr.fit(CH.Southern, mask = CLT_S_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex.S <- secr.fit(CH.Southern, mask = CLT_S_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

save.image("CLTResults2.RData")

#seperately using MS obj
#note need to specify pmix ~ session to est pmix per area (b)
fit.MS.null <- secr.fit(MS, mask = CLT_Both_msk, hcov = "Sex", model = list(D ~ session, lambda0 ~ session, sigma ~ session, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.sexa <- secr.fit(MS, mask = CLT_Both_msk, hcov = "Sex", model = list(D ~ session, lambda0 ~ session+h2, sigma ~ session+h2, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.sexb <- secr.fit(MS, mask = CLT_Both_msk, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)

###################################################################
#refitting to one occasion
#both sessions together

fit.null.both.red <- secr.fit(CH.red, mask = CLT_Both_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.sex.both.red <- secr.fit(CH.red, mask = CLT_Both_msk, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

save.image("CLTResultsRed.RData")

#seperately using MS obj
#note need to specify pmix ~ session to est pmix per area (b)
fit.MS.null.red <- secr.fit(MS.red, mask = CLT_Both_msk, hcov = "Sex", model = list(D ~ session, lambda0 ~ session, sigma ~ session, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.sexa.red <- secr.fit(MS.red, mask = CLT_Both_msk, hcov = "Sex", model = list(D ~ session, lambda0 ~ session+h2, sigma ~ session+h2, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.sexb.red <- secr.fit(MS.red, mask = CLT_Both_msk, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)

save.image("CLTResultsRed.RData")
