library(secrdesign)

# Simple example for secrdesign's GAoptim function (from ?GAoptim)

msk <- make.mask(type = 'rectangular', spacing = 10, nx = 30, ny = 20, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(10,10), spacing = 10)

set.seed(123)

opt <- GAoptim(msk, alltrps, ntraps = 20, detectpar = list(lambda0 = 0.5, sigma = 20), 
               detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, verbose = 1)

plot(msk)
plot(alltrps, detpar = list(pch = 1, cex = 0.4, fg = "gray"), add = TRUE)
plot(opt$optimaltraps, add = TRUE)
minnrRSE(opt, distribution = 'binomial')

# Now update GAoptim to allow fixed detector locations

source("2StepDesign/GAoptim.R")

## Show get same results if no fixed detectors

set.seed(123)

opt2 <- GAoptim(msk, alltrps, ntraps = 20, detectpar = list(lambda0 = 0.5, sigma = 20), 
                detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, verbose = 1)

plot(msk)
plot(alltrps, detpar = list(pch = 1, cex = 0.4, fg = "gray"), add = TRUE)
plot(opt2$optimaltraps, add = TRUE)
minnrRSE(opt2, distribution = 'binomial')


## Example with fixed detectors as a regular grid

set.seed(123)

# Fixed detectors
## origin should be a integer X such that = (alltraps origin) + X * (alltraps spacing)
## spacing should be a multiple of alltraps spacing (make this a 4-5 * sigma for 
## effective space filling)
fixedtrps <- make.grid(nx = 4, ny = 3, origin = c(50,50), spacing = 50)

plot(msk)
plot(alltrps, detpar = list(pch = 1, cex = 0.4, fg = "gray"), add = TRUE)
plot(fixedtrps, detpar = list(fg = "blue"), add = TRUE)

# Remove fixed detectors from alltraps to give remaining detectors available for selection
keep <- rep(TRUE, nrow(alltrps))
for (i in 1:nrow(fixedtrps)) {
  matches <- alltrps[, 1] == fixedtrps[i, 1] & alltrps[, 2] == fixedtrps[i, 2]
  keep[matches] <- FALSE
}
remtrps <- alltrps[keep, , drop = FALSE]

# Note extra argument fixedtraps, and alltraps = remtrps
opt3 <- GAoptim(msk, alltraps = remtrps, fixedtraps = fixedtrps, ntraps = 20, detectpar = list(lambda0 = 0.5, sigma = 20), 
                detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, verbose = 1)

plot(msk)
plot(alltrps, detpar = list(pch = 1, cex = 0.4, fg = "gray"), add = TRUE)
plot(opt3$fixedtraps, detpar = list(fg = "blue"), add = TRUE)
plot(opt3$selectedtraps, add = TRUE)
plot(opt3$optimaltraps, detpar = list(pch = 21, cex = 1.2), add = TRUE)
minnrRSE(opt3, distribution = 'binomial')

## Example with fixed detectors chosen using a space-filling rule

set.seed(123)

# Use all trap locations as the sf sample frame for a spatially balanced sample
nfixed <- 9
alltrps_sf <- sf::st_as_sf(data.frame(id = seq_len(nrow(alltrps)),
                                      x = alltrps$x,
                                      y = alltrps$y),
                           coords = c("x", "y"))
alltrps_sf_poly <- sf::st_convex_hull(sf::st_union(alltrps_sf))
alltrps_poly_sf <- sf::st_sf(id = 1, geometry = alltrps_sf_poly)
grts_out <- spsurvey::grts(alltrps_poly_sf, n_base = nfixed, projcrs_check = FALSE)
fixed_idx <- sf::st_nearest_feature(grts_out$sites_base, alltrps_sf)
fixedtrps_sf <- alltrps[fixed_idx, ]
remtrps_sf <- alltrps[-fixed_idx, ]

plot(msk)
plot(alltrps, detpar = list(pch = 1, cex = 0.4, fg = "gray"), add = TRUE)
plot(fixedtrps_sf, detpar = list(fg = "forestgreen"), add = TRUE)

opt4 <- GAoptim(msk, alltraps = remtrps_sf, fixedtraps = fixedtrps_sf, ntraps = 20,
                detectpar = list(lambda0 = 0.5, sigma = 20),
                detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, verbose = 1)

plot(msk)
plot(alltrps, detpar = list(pch = 1, cex = 0.4, fg = "gray"), add = TRUE)
plot(opt4$fixedtraps, detpar = list(fg = "forestgreen"), add = TRUE)
plot(opt4$selectedtraps, add = TRUE)
plot(opt4$optimaltraps, detpar = list(pch = 21, cex = 1.2), add = TRUE)
minnrRSE(opt4, distribution = 'binomial')
