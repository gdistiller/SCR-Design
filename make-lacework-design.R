library(dplyr)
library(raster)
library(secr)
library(sf)
library(ggplot2)
library(secrdesign)

# function to subset lacework designs to contain N points
source("crop_to_n.R")

# set spacings
smallspacing <- 600
bigspacing <- 2400

# need to set a region in which the design is generated
# region can be a bounding box of mask or a polygon (see ?make.lacework)
temptrap <- make.grid(nx = 6, ny = 6, spacing = 3000)
tempmask <- make.mask(temptrap, buffer = 0, spacing = 2000)
region = data.frame(x = c(min(tempmask$x), min(tempmask$x), max(tempmask$x), max(tempmask$x)),
                    y = c(min(tempmask$y), max(tempmask$y), max(tempmask$y), min(tempmask$y)))

# make the design 
set.seed(4321)
lwrot <- 45
lw <- make.lacework(region = region, 
                    spacing = c(bigspacing, smallspacing),  
                    rotate = lwrot, 
                    detector = "count", 
                    keep.design = TRUE)

plot(tempmask)
plot(lw, add = TRUE)

## Manually remove some traps so end up with the same number of traps as other designs
lw <- crop_to_n(lw$x, lw$y, N = 40, xy_ratio = 2)
lw <- read.traps(data = lw$points, detector = "count")

plot(tempmask)
plot(lw, add = TRUE)
