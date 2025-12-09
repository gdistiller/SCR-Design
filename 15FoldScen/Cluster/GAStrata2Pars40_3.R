#Oct 2025, code to generate GA designs with 40 traps on 45km regular area
#requires the following objects: mask, traplocs
#using foreach and doMPI to get 500 designs
#Using strata 2 values (large sigma, smaller D and L0)
#includes higher ngen
#due to timeout I renamed the results obj to be GA4 or GA5

rm(list=ls())

#################### INITIALIZE MPI ENVIRONMENT ####################
library(secr)
library(secrdesign)
library(kofnGA)

library(doMPI)
library(foreach)
cl <- startMPIcluster()
registerDoMPI(cl)

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

############Load up required functions###########################################
Proposed.traps <- function(poly, alltraps = NULL, nT, D, sigma, lambda0, mask.buff = NULL, sigma.buff = NULL, grid.spacing = NULL, criterion = 4, n.reps = 1, grid = FALSE, GA.ngen = 20, maxim = FALSE, seed = NULL){
  
  Proposed.Traps <- vector(mode = "list", length = n.reps)
  
  for (i in 1:n.reps){
    if (grid == FALSE){
      if (is.null(mask.buff)) {
        mask.buff <- make.mask(traps = alltraps, type = "trapbuffer", spacing = 2/3 * sigma.buff, buffer = 3*sigma.buff) #creates a buffered mask with spacing 2/3 sigma
      } 
      
      opt <- GAoptimTest(mask = mask.buff, alltraps = alltraps, ntraps = nT, detectpar = list(lambda0 = lambda0, sigma = sigma), 
                         criterion = criterion, detectfn = 'HHN', D = D, noccasions = 1, ngen = GA.ngen, verbose = 0, max = maxim, seed = seed)
      Proposed.Traps[[i]] <- opt$optimaltraps
    } else {
      opt.grid <- grid.design(poly, grid.size = grid.spacing, nT)
      Proposed.Traps[[i]] <- opt.grid
    }
  }
  return(Proposed.Traps)
}

Enr.2pops <- function (v, alltraps, mask, detectpar, detectfn, noccasions, 
                       detector, D, crit, weights = c(0.5,0.5), max = FALSE) {
  
  traps <- subset(alltraps, v)
  
  strata = FALSE
  if (length(detectpar[[2]])>1) strata = TRUE
  
  if (is.function(crit)) {
    -crit(D = D, traps = traps, mask = mask, noccasions = noccasions, 
          detectpar = detectpar, detectfn = detectfn)[1]
  }
  if (strata==F){
    if (crit<5) {
      enrm <- Enrm(D = D, traps = traps, mask = mask, noccasions = noccasions, 
                   detectpar = detectpar, detectfn = detectfn)
      c(-enrm[1], -enrm[2], -enrm[3], -(min(enrm[1],enrm[2])))[crit]    
    }
    else {
      en2 <- En2(D = D, traps = traps, mask = mask, noccasions = noccasions, 
                 detectpar = detectpar, detectfn = detectfn)
      c(-en2[2], -sum(en2))[crit-4]
    } 
  } else {  #for 2 populations
    if (crit<5) {
      enrm1 <- Enrm(D = D[1], traps = traps, mask = mask, noccasions = noccasions, 
                    detectpar = list(lambda0 = detectpar[[1]][1], sigma = detectpar[[2]][1]), detectfn = detectfn)
      enrm2 <- Enrm(D = D[2], traps = traps, mask = mask, noccasions = noccasions, 
                    detectpar = list(lambda0 = detectpar[[1]][2], sigma = detectpar[[2]][2]), detectfn = detectfn)
      
      if (max==F){
        obj.n <- -(weights[1]*enrm1[1] + weights[2]*enrm2[1])
        obj.r <- -(weights[1]*enrm1[2] + weights[2]*enrm2[2])
        obj.m <- -(weights[1]*enrm1[3] + weights[2]*enrm2[3])
        obj.min <- -(weights[1]*min(enrm1[1],enrm1[2]) + weights[2]*min(enrm2[1],enrm2[2]))
        
        c(obj.n, obj.r, obj.m, obj.min)[crit]  
      } else {
        obj.n <- -(min(enrm1[1], enrm2[1]))
        obj.r <- -(min(enrm1[2], enrm2[2]))
        obj.m <- -(min(enrm1[3], enrm2[3]))
        obj.min <- -(min(min(enrm1[1],enrm1[2]), min(enrm2[1],enrm2[2])))
        
        c(obj.n, obj.r, obj.m, obj.min)[crit]  
      }
    }
    else {
      en21 <- En2(D = D[1], traps = traps, mask = mask, noccasions = noccasions, 
                  detectpar = list(lambda0 = detectpar[[1]][1], sigma = detectpar[[2]][1]), detectfn = detectfn)
      en22 <- En2(D = D[2], traps = traps, mask = mask, noccasions = noccasions, 
                  detectpar = list(lambda0 = detectpar[[1]][2], sigma = detectpar[[2]][2]), detectfn = detectfn)
      
      if (max==F){
        obj.n2 <- -(weights[1]*en21[2]+weights[2]*en22[2])
        obj.pb <- -(weights[1]*sum(en21)+weights[2]*sum(en22))
        
        c(obj.n2, obj.pb)[crit-4] 
      } else {
        obj.n2 <- -(min(en21[2], en22[2]))
        obj.pb <- -(min(sum(en21), sum(en22)))
        
        c(obj.n2, obj.pb)[crit-4]
      }
    }
  }
}

GAoptimTest <- function(mask, alltraps, ntraps, detectpar, noccasions, detectfn = c("HHN", "HHR", "HEX", "HAN", "HCG"),
                        D = NULL, criterion = 4, weights = c(0.5, 0.5), seed = NULL, ...){
  
  detectfn <- match.arg(detectfn)
  
  ## criterion (1 = En, 2 = Er, 3 = Em, 4 = min(En,Er), 5 = En2, 6 = En+En2)
  if (!is.function(criterion) && (criterion<1 || criterion>6)) stop ("invalid criterion code")
  
  if(missing(mask)) stop("Must supply a 'mask' object (coords of the study area)")
  if(missing(alltraps))   stop("Must supply a 'traps' object (all possible trap locations)")
  
  if (!inherits(mask, "mask")) stop ("mask should be a mask object")
  if (!inherits(alltraps, "traps")) stop ("alltraps should be a traps object")    
  
  detector <- match.arg(detector(alltraps), choices = c("count", "proximity", "multi"))
  if (noccasions == 1 && detector == "multi") stop ("multi detector requires > 1 occasion")
  if(!is.null(seed)) set.seed(seed)
  
  if (ms(mask) || ms(traps)) stop ("mask and traps should be single-session")
  
  #---------------------------------------------------------------------------
  
  des <- kofnGA::kofnGA(n = nrow(alltraps), 
                        k  = ntraps, 
                        OF = Enr.2pops,
                        ...,
                        alltraps    = alltraps,
                        mask        = mask,
                        detectpar   = detectpar,
                        noccasions  = noccasions,
                        detectfn    = detectfn,
                        detector    = detector,
                        D           = if (is.null(D)) 1 else D,
                        crit        = criterion)
  
  optimaltraps <- subset(alltraps, des$bestsol)
  
  if (!is.null(D)) {
    if (length(D)>1){
      optimalenrm1 <- Enrm(D = D[1], traps = optimaltraps, mask = mask, 
                           noccasions = noccasions, detectpar = list("lambda0" = detectpar[[1]][1], "sigma" = detectpar[[2]][1]), detectfn = detectfn)
      optimalenrm2 <- Enrm(D = D[2], traps = optimaltraps, mask = mask, 
                           noccasions = noccasions, detectpar = list("lambda0" = detectpar[[1]][1], "sigma" = detectpar[[2]][2]), detectfn = detectfn)
      optimalenrm <- c(optimalenrm1,optimalenrm2)
    } else {
      optimalenrm <- Enrm(D = D, traps = optimaltraps, mask = mask, 
                          noccasions = noccasions, detectpar = detectpar, detectfn = detectfn)
    }
  }
  else {
    optimalenrm <- NULL
  }
  
  out <- list(
    mask         = mask, 
    alltraps     = alltraps, 
    detectpar    = detectpar, 
    noccasions   = noccasions,
    detectfn     = detectfn,
    D            = D,
    criterion    = criterion,
    des          = des, 
    optimaltraps = optimaltraps,
    optimalenrm  = optimalenrm
    ## do not include minnrRSE - it depends on extra arguments CF, distribution
  )
  
  class(out) <- "GAoptim"
  out
  
}
###################################################################
load("SCRObjs.RData")

mask <- res.objs[[1]] ; trap.locs <- res.objs[[2]]

#Set values for both strata
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor ; L0 <- c(L01,L02)
D1 <- 0.05 ; D2 <- D1 / DiffFactor ; D <- c(D1, D2)
nT <- 40 ; nreps = 500

###################################
#get proposed designs using strata 1
###################################

#GA4.Strata2Pars3.40.list <- vector("list", nreps)
GA5.Strata2Pars3.40.list <- vector("list", nreps)

#using foreach 
#GA4.Strata2Pars3.40.list <- foreach (r = 1:nreps, .combine='c') %dopar% {Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D2, sigma = sigma2, 
#                                                                  lambda0 = L02, sigma.buff = sigma2, grid.spacing = NULL, criterion = 4, 
#                                                                  n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 500)}
#save.image("GAStrata2Pars3.RData")
#renamed to be GA4Strata2Pars3.RData

GA5.Strata2Pars3.40.list <- foreach (r = 1:nreps, .combine='c') %dopar% {Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D2, sigma = sigma2, 
                                                                  lambda0 = L02, sigma.buff = sigma2, grid.spacing = NULL, criterion = 5, 
                                                                  n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 500)}
save.image("GA5Strata2Pars3.RData")

###################################################################

closeCluster(cl)
mpi.quit(save = "no")