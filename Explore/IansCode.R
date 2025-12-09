#Sept 2023, bits of Ian's code

############################
#generating optimal designs
############################
#use GAoptim to get proposed designs
#regular 2 sigma and optimised grid follow

# user parameters
lambda0 <- 2 
dens_per_100km2 <- 2 # mean animal density per 100km2, SLs are ~1
D <- dens_per_100km2 / 10000
sigma <- 3000
b_ac <- 1

nT <- 20
buffer <- 4 * sigma
dt <- "count"

# an artificial example, note the bottom-left corner is at origin
msk <- make.mask(type = 'rectangular', spacing = 1000, nx = 30, ny = 20, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(1000,1000), spacing = 1000, detector = "count")

# 50 generations for demonstration, use more in practice
opt <- GAoptim(mask = msk, alltraps = alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
               detectfn = 'HHN', D = D, noccasions = 1, ngen = 5, verbose = 1)

plot(msk)
plot(opt$optimaltraps, add = TRUE)
minnrRSE(opt, distribution = 'binomial')

# Using a criterion function
# En2 is unsuitable as a criterion function as it returns 2 values
# This function selects the second as the (unique) criterion
fn <- function(...) En2(...)[2]
opt2 <- GAoptim(msk, alltrps, ntraps = 20, detectpar = list(lambda0 = lambda0, sigma = sigma), 
                detectfn = 'HHN', D = D, noccasions = 1, ngen = 3, verbose = 1, criterion = fn)

plot(msk)
plot(opt2$optimaltraps, add = TRUE)
minnrRSE(opt2, distribution = 'binomial')

# grid design under constant D
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(stringr)
library(purrr)
library(secr)
library(secrdesign)
library(oSCR)
library(raster)
library(gdistance)
library(kofnGA)

# create mask and then a buffered mask (the mask used for the designs)
mask <- secr::raster(msk)
mask_df <- data.frame(coordinates(mask))
mask <- read.mask(data = mask_df)
plot(mask, dots=F)

buffer_mult <- ceiling(buffer / cellsize)
newmask_df <- expand.grid(x = seq(from = min(mask_df$x) - buffer_mult * cellsize, 
                                  to = max(mask_df$x) + buffer_mult * cellsize, 
                                  by = cellsize),
                          y = seq(from = min(mask_df$y) - buffer_mult * cellsize, 
                                  to = max(mask_df$y) + buffer_mult * cellsize, 
                                  by = cellsize)) 
newmask_df <- newmask_df %>% 
  mutate(keep = (apply(e2dist(as.matrix(newmask_df), as.matrix(mask_df)), 1, min) <= buffer)) %>%
  filter(keep == TRUE) %>% dplyr::select(-keep)
newmask <- read.mask(data = newmask_df)
plot(newmask,axes=T)

# create grid of possible trap locations
# reduce resolution of mesh so have fewer possible camera locations
cellsize <- 2/3 * sigma
red_factor <- cellsize[1] / attr(msk, "spacing")
if ((trunc(red_factor) - red_factor) != 0) stop("Check spacing, causing non-integer reduction factor for mesh")

alltraps_df <- data.frame(coordinates(data.frame(alltrps)))
alltraps <- as.matrix(alltraps_df)
plot(read.mask(data = alltraps_df), col = 'blue', add = TRUE, axes =T)  

#####################################
grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                             x = as.numeric(), y = as.numeric(), trap_id = as.integer())
opt_grid_traps_all <- data.frame(nT = as.integer(), sigma = as.numeric(), beta0 = as.numeric(), buffer = as.numeric(), dt = as.character(),
                                 x = as.numeric(), y = as.numeric(), trap_id = as.integer())
######################
#regular 2 sigma grid
######################

for(n in c(5,10,15)){
  
  # hacky way to make a polygon just bigger than boundary of mask points
  # used later to decide which randomly generated detectors are on the mask and so allowed
  sap_1 <- mask_df %>% st_as_sf(coords = c("x", "y")) %>% st_buffer(dist = cellsize, endCapStyle = "SQUARE") %>% st_union() 
  sap_1 <- sap_1 %>% st_buffer(dist = -cellsize * 0.4, endCapStyle = "SQUARE")
  
  # place a grid over the area, with cells 2 * sigma apart
  my_grid <- st_make_grid(sap_1, cellsize = c(2 * sigma, 2 * sigma), what = "centers") %>%
    st_intersection(sap_1)
  grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
  
  # choose a subset of nT of these, starting from a random trap and then choosing nT nearest neighbours
  all_grid_traps <- list()
  
  grid_traps <- grid_traps_full
  set.seed(700)
  xy_rand <- grid_traps[sample(1:nrow(grid_traps), 1), ]
  grid_traps$dist2pt <- (grid_traps$x - xy_rand$x)^2 + (grid_traps$y - xy_rand$y)^2
  grid_traps <- grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
  grid_traps_all <- rbind(grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, grid_traps))
  
  ######################
  #optimised grid
  ######################
  grid_traps <- read.traps(data = grid_traps, detector = dt)
  
  # compute optimal spacing for previous grid 
  dd <- optimalSpacing(D = D, traps = grid_traps, detectpar = list(lambda0 = lambda0, sigma = sigma), noccasions = 1)$rotRSE$optimum.spacing
  
  # sometimes optimal spacing is too big to fit desired number of traps in; if so, reduce dd until it is
  n_opt_grid <- 0
  red_dd <- 0.9
  
  while(n_opt_grid < n){
    
    dd <- dd * red_dd
    
    # place a grid over the area, with cells optimally spaced
    my_grid <- st_make_grid(sap_1, cellsize = c(dd, dd), what = "centers") %>% st_intersection(sap_1)
    opt_grid_traps_full <- st_coordinates(my_grid) %>% as.data.frame() %>% rename(x = X, y = Y)
    
    # choose a subset of nT of these points, starting from a random trap and then choosing nT nearest neighbours
    opt_grid_traps <- opt_grid_traps_full
    set.seed(700)
    xy_rand <- opt_grid_traps[sample(1:nrow(opt_grid_traps), 1), ]
    opt_grid_traps$dist2pt <- (opt_grid_traps$x - xy_rand$x)^2 + (opt_grid_traps$y - xy_rand$y)^2
    opt_grid_traps <- opt_grid_traps %>% filter(rank(dist2pt, ties.method = "random") <= n) %>% mutate(trap_id = 1)
    
    n_opt_grid <- nrow(opt_grid_traps)
    
  }
  
  opt_grid_traps_all <- rbind(opt_grid_traps_all, cbind(nT = n, lambda0 = lambda0, sigma = sigma, buffer = buffer, dt = dt, opt_grid_traps))
  
}
