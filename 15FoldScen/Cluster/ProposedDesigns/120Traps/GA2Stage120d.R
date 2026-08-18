#June 26, code to generate 2 stage designs with 120 traps on 45km regular area
#version d uses En2 for 2nd stage
#requires the following objects: mask, traplocs
#using foreach and doMPI to get 500 designs

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

# default pen_fn
GApenfn <- function (traps, sigma) {
  # find out how many detector pairs are between 2.5-3.5 and 3.5-4.5 sigma apart
  breaks <- c(0, 2.499, 3.499, 4.499, Inf) * sigma     # why 4.449? assume typo
  d <- as.matrix(dist(traps))  # for compatibility
  tabulate(cut(d, breaks = breaks))[2:3]
}

# used within new fn
compactSample <- function (traps, n) {
  # closest ntraps to random start
  xy_rand <- traps[sample.int(nrow(traps), 1), ]  
  dist2pt <- distancetotrap(traps, xy_rand)
  OK <- rank(dist2pt, ties.method = "random") <= n
  subset(traps, OK)
}

# used within new fn
as_traps_object <- function (x, detector) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "traps")) {
    if (detector(x) != detector) stop("'fixedtraps' detector type should match 'alltraps'")
    return(x)
  }
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("'fixedtraps' should be NULL, a matrix/data.frame of coordinates, or a traps object")
  }
  if (ncol(x) < 2) stop("'fixedtraps' should have at least two columns (x, y)")
  x <- as.data.frame(x[, 1:2, drop = FALSE])
  names(x) <- c("x", "y")
  read.traps(data = x, detector = detector)
}

# used within new fn
combine_traps <- function (selectedtraps, fixedtraps) {
  if (is.null(fixedtraps) || nrow(fixedtraps) == 0) selectedtraps
  else if (nrow(selectedtraps) == 0) fixedtraps
  else rbind(selectedtraps, fixedtraps)
}

# Objective function
OF <- function (v, alltraps, fixedtraps, mask, detectpar, detectfn, noccasions, detector, D, crit, penalty, g_penvector) {
  
  selectedtraps <- subset(alltraps, v)
  traps <- combine_traps(selectedtraps, fixedtraps)
  
  # penalty for too clustered
  if (!is.null(penalty)) {
    penvector <- penalty$pen_fn(traps, detectpar$sigma)
    penalty <- penalty$pen_wt * sum(pmax(0, g_penvector - penvector))
  }
  else {
    penalty <- 0
  } 
  
  if (length(detectpar$lambda0) > 1)
    stop ("this implementation does not allow varying lambda0")
  
  if (is.function(crit)) {
    -crit(D = D, traps = traps, mask = mask, noccasions = noccasions, 
          detectpar = detectpar, detectfn = detectfn)[1]
  }
  else if (crit<5) {
    enrm <- Enrm(D = D, traps = traps, mask = mask, noccasions = noccasions, 
                 detectpar = detectpar, detectfn = detectfn)
    c(-enrm[1], -enrm[2], -enrm[3], penalty-(min(enrm[1],enrm[2])))[crit]    
  }
  else {
    en2 <- En2(D = D, traps = traps, mask = mask, noccasions = noccasions, 
               detectpar = detectpar, detectfn = detectfn)
    -c(en2[2], sum(en2))[crit-4]
  }
}

#extended GAoptim
GAoptim <- function(mask, alltraps, ntraps, fixedtraps = NULL, detectpar, noccasions, detectfn  = c("HHN", "HHR", "HEX", "HAN", "HCG"), D = NULL, criterion = 5, penalty = NULL, seed = NULL, ...){
  
  detectfn <- match.arg(detectfn)
  
  ## criterion (1 = En, 2 = Er, 3 = Em, 4 = min(En,Er), 5 = En2, 6 = En+En2)
  if (!is.function(criterion) && (criterion<1 || criterion>6)) stop ("invalid criterion code")
  
  if(missing(mask)) stop("Must supply a 'mask' object (coords of the study area)")
  if(missing(alltraps))   stop("Must supply a 'traps' object (all possible trap locations)")
  
  if (!inherits(mask, "mask")) stop ("mask should be a mask object")
  if (!inherits(alltraps, "traps")) stop ("alltraps should be a traps object")    
  
  detector <- match.arg(detector(alltraps), choices = c("count", "proximity", "multi"))
  fixedtraps <- as_traps_object(fixedtraps, detector = detector)
  nfixedtraps <- if (is.null(fixedtraps)) 0 else nrow(fixedtraps)
  nselect <- ntraps - nfixedtraps
  if (nselect < 0) stop("'ntraps' must be at least nrow(fixedtraps)")
  if (nselect > nrow(alltraps)) stop("'ntraps - nrow(fixedtraps)' exceeds nrow(alltraps)")
  if (noccasions == 1 && detector == "multi") stop ("multi detector requires > 1 occasion")
  if(!is.null(seed)) set.seed(seed)
  
  if (ms(mask) || ms(alltraps) || (!is.null(fixedtraps) && ms(fixedtraps))) {
    stop ("mask and traps should be single-session")
  }
  
  #---------------------------------------------------------------------------
  if (!is.null(penalty)) {
    # penalty reference vector (Durbach et al. 2021)
    # use default penalty function (see above) if none provided
    if (is.null(penalty$pen_fn)) penalty$pen_fn <- GApenfn  
    
    # find distribution of trap spacings on a close to regular grid, to ensure 
    # later optimized grid has spaced enough detectors sufficiently far apart 
    # to get low var(sigma)
    
    # polygon to represent region of interest
    pg <- st_union(gridCells(alltraps))
    # place a grid over the area, with cells pen_gridsigma * sigma apart
    cellsize <- penalty$pen_gridsigma * detectpar$sigma
    grid_traps <- make.systematic(region = pg, spacing = cellsize)
    grid_traps <- compactSample(grid_traps, nselect)
    grid_traps <- combine_traps(grid_traps, fixedtraps)
    # target vector (e.g., number of traps in each distance bracket)
    g_penvector <- penalty$pen_fn(grid_traps, detectpar$sigma)
  }
  else {
    g_penvector <- NA 
  }
  #---------------------------------------------------------------------------
  
  if (nselect == 0) {
    des <- list(bestsol = rep(FALSE, nrow(alltraps)))
  }
  else {
    des <- kofnGA::kofnGA(n = nrow(alltraps), 
                          k  = nselect, 
                          OF = OF,
                          ...,
                          alltraps    = alltraps,
                          fixedtraps  = fixedtraps,
                          mask        = mask,
                          detectpar   = detectpar,
                          noccasions  = noccasions,
                          detectfn    = detectfn,
                          detector    = detector,
                          D           = if (is.null(D)) 1 else D,
                          crit        = criterion,
                          penalty     = penalty,
                          g_penvector = g_penvector
    )
  }
  
  selectedtraps <- subset(alltraps, des$bestsol)
  optimaltraps <- combine_traps(selectedtraps, fixedtraps)
  
  if (!is.null(D)) {
    optimalenrm <- Enrm(D = D, traps = optimaltraps, mask = mask, 
                        noccasions = noccasions, detectpar = detectpar, detectfn = detectfn)
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
    penalty      = penalty,
    criterion    = criterion,
    fixedtraps   = fixedtraps,
    des          = des, 
    selectedtraps = selectedtraps,
    optimaltraps = optimaltraps,
    optimalenrm  = optimalenrm
    ## do not include minnrRSE - it depends on extra arguments CF, distribution
  )
  
  class(out) <- "GAoptim"
  out
}

#function to generate designs following the 2 stage approach
#a space filling approach taken for the first half 
#optimise 2nd stage based on S1 values (using En2)
#all traps must be all possible trap locations as a traps / df obj
Two.stage.design <- function(msk, alltraps, L0, Sig, Dens, ntrps, nfixed, ngens, nreps = 1){
  
  Proposed.Traps <- vector(mode = "list", length = nreps)
  Stage1.Traps <- vector(mode = "list", length = nreps)
  Stage2.Traps <- vector(mode = "list", length = nreps)
  
  #construct objs for all possible trap locations
  alltrps_sf <- sf::st_as_sf(data.frame(id = seq_len(nrow(alltraps)),
                                        x = alltraps$x,
                                        y = alltraps$y),
                             coords = c("x", "y"))
  alltrps_sf_poly <- sf::st_convex_hull(sf::st_union(alltrps_sf))
  alltrps_poly_sf <- sf::st_sf(id = 1, geometry = alltrps_sf_poly)
  
  #repeat for each design
  for (i in 1:nreps){
    grts_out <- spsurvey::grts(alltrps_poly_sf, n_base = nfixed, projcrs_check = FALSE)
    fixed_idx <- sf::st_nearest_feature(grts_out$sites_base, alltrps_sf)
    fixedtrps_sf <- alltraps[fixed_idx, ]
    remtrps_sf <- alltraps[-fixed_idx, ]
    
    opt <- GAoptim(msk, alltraps = remtrps_sf, fixedtraps = fixedtrps_sf, ntraps = ntrps,
                   detectpar = list(lambda0 = L0, sigma = Sig),
                   detectfn = 'HHN', D = Dens, noccasions = 1, ngen = ngens, verbose = 0)
    
    Proposed.Traps[[i]] <- opt$optimaltraps
    Stage1.Traps[[i]] <- opt$fixedtraps
    Stage2.Traps[[i]] <- opt$selectedtraps
  }
  res <- list("Proposed design" = Proposed.Traps, "Stage1 traps" = Stage1.Traps, "Stage2 traps" = Stage2.Traps)
  return(res)
}

###################################################################
load("SCRObjs.RData")

mask <- SCR.objs$Full$`SCR Mask` ; trap.locs <- SCR.objs$Full$`SCR trap locs`

#Set values for both strata
DiffFactor <- 15
sigma1 = 200 ; sigma2 = sigma1*DiffFactor ; sigma <- c(sigma1, sigma2)
L01 <- 2 ; L02 <- L01/DiffFactor ; L0 <- c(L01, L02)
D1 <- 0.05 ; D2 <- D1 / DiffFactor ; D <- c(D1, D2)
nT <- 120 ; nreps = 500

#####################
#get proposed designs
#####################
TwoStage.120.list <- vector("list", nreps)

TwoStage.120.list <- foreach (r = 1:nreps, .combine='c') %dopar% {Two.stage.design(msk = mask, alltraps = trap.locs, L0 = L01, Sig = sigma1, Dens = D1, ntrps = nT, nfixed = nT/2, ngens = 500)}

save.image("GA2Stage120d.RData")

###################################################################

closeCluster(cl)
mpi.quit(save = "no")
