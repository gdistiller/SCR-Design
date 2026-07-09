# example lacework design

library(secr)

source("make_lacework_fixed_ndet.R")
environment(make.lacework) <- asNamespace("secr")

# set spacings
# unclear what these should be, but:
# smallspacing should probably be smaller than cluster design within-cluster spacing, because
# lacework doesn't give the same regular grid arrangement.
# similarly bigspacing should probably be bigger than cluster design between-cluster spacing, because
# lacework gives plenty of opportunities for recaps at intermediate distances.

smallspacing <- 400
bigspacing <- 6000

# need to set a region in which the design is generated
# region can be a bounding box of mask or a polygon (see ?make.lacework)
temptrap <- make.grid(nx = 20, ny = 20, spacing = 2000)
tempmask <- make.mask(temptrap, buffer = 0, spacing = 2000)
region = data.frame(x = c(min(tempmask$x), min(tempmask$x), max(tempmask$x), max(tempmask$x)),
                    y = c(min(tempmask$y), max(tempmask$y), max(tempmask$y), min(tempmask$y)))

# unconstrained, just removed ndet < 5000 restriction
set.seed(4321)
lwrot <- 45
lw <- make.lacework(region = region, 
                    spacing = c(bigspacing, smallspacing),  
                    rotate = lwrot, 
                    detector = "count", 
                    radius = NULL,
                    keep.design = TRUE)

plot(tempmask)
plot(lw, add = TRUE)

# restrict the number of intersections (with idea these will/may become clusters)
set.seed(4321)
lwrot <- 45
lw <- make.lacework(region = region, 
                    spacing = c(bigspacing, smallspacing),  
                    rotate = lwrot, 
                    nintersections = 6,
                    detector = "count", 
                    radius = NULL,
                    keep.design = TRUE)

plot(tempmask)
plot(lw, add = TRUE)

# further restrict the number of detectors by varying radius internally
set.seed(4321)
lwrot <- 45
lw <- make.lacework(region = region, 
                    spacing = c(bigspacing, smallspacing),  
                    rotate = lwrot, 
                    nintersections = 10,
                    ndetectors = 40,
                    detector = "count", 
                    radius = NULL,
                    keep.design = TRUE)

plot(tempmask)
plot(lw, add = TRUE)

