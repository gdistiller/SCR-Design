#August 2024, revised set of functions used in SCRDesign project
#March 26, added optimal spacing

###################################################################################
#generates the proposed traps without doing Enrm summaries
#alltraps must be a traps obj
#Note that: mesh of traplocs only required for GA and area polygon for
#a) for irregular areas to create the mask and b) for grids   
#buffered mask required for GA optimisation, created from traplocs if not supplied
#set buffer sigma and grid cell resolution seperately
#GAoptim fn will handle SCR pars for two sexes
Proposed.traps <- function(poly, alltraps = NULL, nT, D, sigma, lambda0, mask.buff = NULL, sigma.buff = NULL, grid.spacing = NULL, criterion = 4, n.reps = 1, grid = FALSE, GA.ngen = 20, maxim = FALSE){
  
  Proposed.Traps <- vector(mode = "list", length = n.reps)
  
  for (i in 1:n.reps){
    if (grid == FALSE){
      if (is.null(mask.buff)) {
        mask.buff <- make.mask(traps = alltraps, type = "trapbuffer", spacing = 2/3 * sigma.buff, buffer = 3*sigma.buff) #creates a buffered mask with spacing 2/3 sigma
      } 
      
      opt <- GAoptimTest(mask = mask.buff, alltraps = alltraps, ntraps = nT, detectpar = list(lambda0 = lambda0, sigma = sigma), 
                         criterion = criterion, detectfn = 'HHN', D = D, noccasions = 1, ngen = GA.ngen, verbose = 0, max = maxim)
      Proposed.Traps[[i]] <- opt$optimaltraps
    } else {
      opt.grid <- grid.design(poly, grid.size = grid.spacing, nT)
      Proposed.Traps[[i]] <- opt.grid
    }
  }
  return(Proposed.Traps)
}

#function to generate designs following the 2 stage approach
#a space filling approach taken for the first half based on S1
#optimise 2nd stage based on S2 values (using min(n,r))
#all traps must be all possible trap locations as a traps / df obj
Two.stage.design <- function(msk, alltraps, L0, Sig, Dens, ntrps, nfixed, ngens, nreps){
  
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

#function to create a mask and trap locs objects, includes SCR mask and traps obj and a SF polygon for grids
#sig.factor determines the extent ito sigma, buffer is equal to 3 sigma
#res is the resolution of both
create.extent <- function(sigma, sig.factor = 15, buff.factor = 3, res = 200){
  Extent <- sig.factor * sigma ; ncolrow1 = Extent/res 
  BuffExtent1 <- buff.factor * sigma ; BuffExtent2 <- Extent-BuffExtent1 ; ncolrow2 = (BuffExtent2-BuffExtent1)/res
  
  mask.r <- rast(ncol = ncolrow1, nrow = ncolrow1, xmin = 0, xmax = Extent, ymin = 0, ymax = Extent)
  traplocs.r <- rast(ncol = ncolrow2, nrow = ncolrow2, xmin = BuffExtent1, xmax = BuffExtent2, ymin = BuffExtent1, ymax = BuffExtent2) 
  
  values(mask.r) <- 1:ncell(mask.r) #numbers each cell, seems needed
  values(traplocs.r) <- 1:ncell(traplocs.r) 
  
  #create polygon
  msk.p1 <- rbind(c(0, 0), c(Extent, 0), c(Extent, Extent), c(0, Extent), c(0,0))
  mask.v <- vect(msk.p1, type = "polygons")
  
  traplocs.p1 <- rbind(c(BuffExtent1, BuffExtent1), c(BuffExtent2, BuffExtent1), c(BuffExtent2, BuffExtent2), c(BuffExtent1, BuffExtent2), c(BuffExtent1,BuffExtent1))
  traplocs.v <- vect(traplocs.p1, type = "polygons")
  traplocs.sf <- sf::st_as_sf(traplocs.v)
  
  #create mask and traps obj
  #get coordinates
  msk.coords <- crds(mask.r, df = T)
  traploc.coords <- crds(traplocs.r, df=T)
  
  mask <- read.mask(data = msk.coords, spacing = res)
  trap.locs <- read.traps(data = traploc.coords, detector = "count", spacing = res)
  objs <- list("SCR Mask" = mask, "SCR trap locs" = trap.locs, "SF traps polygon" = traplocs.sf)
}

#new OF to pass to GAOptim
#allows one to weight Enrm for two strata
#this version allows one to maximise the min of the two pops (so focus on worst performing)
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

#version of GAoptim that uses the OF above
#this one uses a common lambda0
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

#function to randomly generate a grid design
#poly has to be an sf or sfc object, shape file is a sf object
#randomly picks a grid point and selects nearest neighbours
grid.design <- function(poly, grid.size, nT){
  
  # place a grid over the area, with cells 2 * sigma apart
  my_grid <- st_make_grid(poly, cellsize = grid.size, what = "centers") %>%
    st_intersection(poly)
  grid_traps <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= nT) %>% mutate(trap_id = 1)
  opt.grid <- read.traps(data = grid_traps[,1:2], detector = "count")
  return(opt.grid)
}

#' Generate detector coordinates for a clustered grid design
#' Uses a random grid of cluster centroids
#' @param centroids A matrix or data frame with columns 'x' and 'y' (or first two columns used)
#' @param spacing Distance between detectors within a cluster
#' @param n_rows Number of rows in each cluster grid
#' @param n_cols Number of columns in each cluster grid
#' @return A data frame with columns: cluster, trap, x, y

make_cluster_traps <- function(centroids, spacing = 100, n_rows = 2, n_cols = 2) {
  
  # Build the local grid offsets, centred on (0, 0)
  col_offsets <- (seq_len(n_cols) - (n_cols + 1) / 2) * spacing
  row_offsets <- (seq_len(n_rows) - (n_rows + 1) / 2) * spacing
  
  local_grid <- expand.grid(x = col_offsets, y = row_offsets)
  
  # Ensure centroids is a data frame with x/y columns
  if (is.matrix(centroids)) {
    centroids <- as.data.frame(centroids)
  }
  colnames(centroids)[1:2] <- c("x", "y")
  
  # Expand each centroid by the local grid offsets
  n_clusters  <- nrow(centroids)
  n_per_cluster <- nrow(local_grid)
  
  result <- do.call(rbind, lapply(seq_len(n_clusters), function(i) {
    data.frame(
      cluster = i,
      trap    = seq_len(n_per_cluster),
      x       = centroids$x[i] + local_grid$x,
      y       = centroids$y[i] + local_grid$y
    )
  }))
  
  rownames(result) <- NULL
  result
}

#sim function
#proposed designs need to be passed as list
#raw = T keeps raw datasets and summary stats of data
#order of groups is F then M
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


#function to extract and summarise estimates for two sexes
#coded for F to come first
summ.results <- function(results.obj, par = "D", true.values, numests, ylim.ceiling, mean = TRUE){
  F.ests <- vector("numeric", numests)
  M.ests <- vector("numeric", numests)
  
  for (i in 1:numests){
    if (par == "D"){
      F.ests[[i]] <- results.obj$output[[i]][[1]][[1]][1,2]
      M.ests[[i]] <- results.obj$output[[i]][[1]][[2]][1,2]
      plot.label1 <- "Female Density"
      plot.label2 <- "Male Density"
    } else {
      if (par=="Sigma"){
        F.ests[[i]] <- results.obj$output[[i]][[1]][[1]][3,2]
        M.ests[[i]] <- results.obj$output[[i]][[1]][[2]][3,2]   
        plot.label1 <- "Female Sigma"
        plot.label2 <- "Male Sigma"
      } else {
        F.ests[[i]] <- results.obj$output[[i]][[1]][[1]][2,2]
        M.ests[[i]] <- results.obj$output[[i]][[1]][[2]][2,2]
        plot.label1 <- "Female Lambda0"
        plot.label2 <- "Male Lambda0"
      }
    }
  }
  
  #calculate RB and RSE
  if (mean == TRUE){
    F.RB <- (mean(F.ests) - true.values[1]) / true.values[1] 
    M.RB <- (mean(M.ests) - true.values[2]) / true.values[2]
  } else {
    F.RB <- (median(F.ests) - true.values[1]) / true.values[1] 
    M.RB <- (median(M.ests) - true.values[2]) / true.values[2]
  }
    
  F.RSE <- sd(F.ests) / mean(F.ests) 
  M.RSE <- sd(M.ests) / mean(M.ests)  
  
  #produce side by side boxplots  
  par(mfrow=c(1,2))
  boxplot(F.ests, main = paste(plot.label1, "estimates", sep =" "), ylim = c(0, ylim.ceiling))
  abline(h = true.values[1], col = 'red')
  
  boxplot(M.ests, main = paste(plot.label2, "estimates", sep =" "), ylim = c(0, ylim.ceiling))
  abline(h = true.values[2], col = 'red')
  res <- list("Rel bias" = c(F.RB, M.RB), "Rel SE" = c(F.RSE, M.RSE), "Female ests" = F.ests, "Male ests" = M.ests)
  print(res[[1]])
  print(res[[2]])
  return(res)
}

##############################################################################
## Two-group optimal spacing (absolute spacing search)
## Poisson only, no simulation block
## Objective: maximise either min(En1 + En2, Er1 + Er2) or min(En1, En2, Er1, Er2)
##############################################################################

optimalSpacing2G <- function(
    D1, D2,
    traps0,
    detectpar1, detectpar2,
    noccasions,
    nrepeats = 1,
    detectfn1 = "HHN",
    detectfn2 = "HHN",
    xsigma = 4,
    spacing_m = seq(200, 4000, 200),
    criterion = c("sum_min", "all_min"),
    CF = 1.0,
    ...
) {
  
  criterion <- match.arg(criterion)
  
  detectfn1 <- secr:::secr_valid.detectfn(detectfn1, valid = c(0, 1, 2, 14:19))
  detectfn2 <- secr:::secr_valid.detectfn(detectfn2, valid = c(0, 1, 2, 14:19))
  dfc1 <- secrdesign:::dfcast(detectfn1, detectpar1)
  dfc2 <- secrdesign:::dfcast(detectfn2, detectpar2)
  detectfn1 <- dfc1$detectfn
  detectfn2 <- dfc2$detectfn
  detectpar1 <- dfc1$detectpar
  detectpar2 <- dfc2$detectpar
  
  base_spacing <- spacing(traps0)
  
  args <- list(...)
  defaultmaskargs <- list(nx = 64, type = "trapbuffer")
  dotsmask <- args[names(args) %in% c("nx", "type", "poly", "poly.habitat")]
  maskargs_base <- secrdesign:::replacedefaults(defaultmaskargs, dotsmask)
  
  if (any(detector(traps0) == "single")) {
    warning("treating single-catch traps as multi-catch", call. = FALSE)
    detector(traps0) <- "multi"
  }
  
  values <- lapply(
    spacing_m,
    getCrit2G_abs,
    traps0,
    base_spacing,
    xsigma,
    detectpar1,
    detectpar2,
    D1,
    D2,
    noccasions,
    nrepeats,
    detectfn1,
    detectfn2,
    maskargs_base,
    CF,
    criterion
  )
  values <- do.call(rbind, values)
  values <- as.data.frame(values)
  names(values) <- c("spacing", "En1", "En2", "Er1", "Er2", "crit")
  
  opt <- if (!is.na(CF) && nrow(values) > 0 && diff(range(values$spacing)) > 0) {
    interpCritMax(values)
  } else {
    list(minimum = NA, objective = NA)
  }
  
  out <- list(
    values = values,
    optimum.spacing = opt$minimum,
    maximum.crit = opt$objective,
    criterion = criterion,
    traps.base = traps0,
    detectpar1 = detectpar1,
    detectpar2 = detectpar2
  )
  class(out) <- c("optimalSpacing2G", "list")
  out
}

##############################################################################

getCrit2G_abs <- function(
    S,
    traps0,
    base_spacing,
    xsigma,
    detectpar1,
    detectpar2,
    D1,
    D2,
    noccasions,
    nrepeats,
    detectfn1,
    detectfn2,
    maskargs_base,
    CF,
    criterion
) {
  
  scalefac <- S / base_spacing
  trapS <- scale_traps(traps0, scalefac)
  
  maskargs <- maskargs_base
  maskargs$traps <- trapS
  maskargs$buffer <- xsigma * max(detectpar1$sigma, detectpar2$sigma) + S
  mask <- do.call(make.mask, maskargs)
  
  En1 <- Enrm(D1, trapS, mask, detectpar1, noccasions, detectfn1) * nrepeats
  En2 <- Enrm(D2, trapS, mask, detectpar2, noccasions, detectfn2) * nrepeats
  Er1 <- En1[2]
  Er2 <- En2[2]
  En1 <- En1[1]
  En2 <- En2[1]
  
  critval <- switch(
    criterion,
    sum_min = min(En1 + En2, Er1 + Er2),
    all_min = min(En1, En2, Er1, Er2)
  )
  
  c(S, En1, En2, Er1, Er2, critval * CF)
}

##############################################################################

interpCritMax <- function(values) {
  f <- approxfun(values$spacing, values$crit, rule = 2)
  opt <- optimize(f, interval = range(values$spacing), maximum = TRUE)
  list(minimum = opt$maximum, objective = opt$objective)
}

##############################################################################

scale_traps <- function(traps, scalefac) {
  tr <- traps
  tr$x <- tr$x * scalefac
  tr$y <- tr$y * scalefac
  for (nm in c("spacing", "spacex", "spacey")) {
    if (!is.null(attr(tr, nm))) {
      attr(tr, nm) <- attr(tr, nm) * scalefac
    }
  }
  tr
}

##############################################################################

