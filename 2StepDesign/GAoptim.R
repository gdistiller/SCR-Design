##############################################################################
## package 'secrdesign'
## GAoptim.R
## adapted from code of Ian Durbach 
## 2022-10-02,06,22
##############################################################################

# default pen_fn
GApenfn <- function (traps, sigma) {
  # find out how many detector pairs are between 2.5-3.5 and 3.5-4.5 sigma apart
  breaks <- c(0, 2.499, 3.499, 4.499, Inf) * sigma     # why 4.449? assume typo
  d <- as.matrix(dist(traps))  # for compatibility
  tabulate(cut(d, breaks = breaks))[2:3]
}

compactSample <- function (traps, n) {
  # closest ntraps to random start
  xy_rand <- traps[sample.int(nrow(traps), 1), ]  
  dist2pt <- distancetotrap(traps, xy_rand)
  OK <- rank(dist2pt, ties.method = "random") <= n
  subset(traps, OK)
}

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

combine_traps <- function (selectedtraps, fixedtraps) {
  if (is.null(fixedtraps) || nrow(fixedtraps) == 0) selectedtraps
  else if (nrow(selectedtraps) == 0) fixedtraps
  else rbind(selectedtraps, fixedtraps)
}

#-------------------------------------------------------------------------------
# Objective function
OF <- function (v, 
                alltraps, 
                fixedtraps,
                mask, 
                detectpar,
                detectfn,
                noccasions, 
                detector, 
                D, 
                crit, 
                penalty,
                g_penvector) {
  
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
#-------------------------------------------------------------------------------

GAoptim <- function(
    mask,
    alltraps, 
    ntraps, 
    fixedtraps = NULL,
    detectpar, 
    noccasions,
    detectfn  = c("HHN", "HHR", "HEX", "HAN", "HCG"),
    D         = NULL,
    criterion = 4,
    penalty   = NULL,
    seed      = NULL,
    ...){
  
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
