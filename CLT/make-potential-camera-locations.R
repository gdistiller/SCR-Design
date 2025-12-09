#Sept 2024, 
#makes potential trap locs and masks for CLT area
#does it for the simnple and complex area (with parts cut out)
#also does it for the two areas seperately

library(sf)
library(leaflet)
library(secr)

x <- readRDS("CLT/clt_surveyregion_anita.Rds")
class(x)

# plot regions to see what's happening
plot(st_geometry(x), axes= TRUE)
plot(st_geometry(x[1,]), axes= TRUE) #simple area without polygons removed
plot(st_geometry(x[2,]), add= TRUE) # Grabouw area to remove
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[3,]), add= TRUE) # small area to remove
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[4,]), add= TRUE) # small area to remove
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[5,]), add= TRUE) # small area to add

# convert to UTM so can take intersections, unions etc
my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"
x1_utm <- x[1,] %>% st_transform(my_crs)
x2_utm <- x[2,] %>% st_transform(my_crs)
x3_utm <- x[3,] %>% st_transform(my_crs)
x4_utm <- x[4,] %>% st_transform(my_crs)
x5_utm <- x[5,] %>% st_transform(my_crs)

# remove areas x2, x3, x4 and add area x5
potential_cams <- x1_utm %>% 
  st_difference(x2_utm) %>% 
  st_difference(x3_utm) %>% 
  st_difference(x4_utm) %>%
  st_union(x5_utm)
plot(potential_cams)

# plot on map
leaflet()%>%
  addProviderTiles(provider= "CartoDB.Positron") %>%
  addPolygons(data = potential_cams %>% st_transform(4326), color = "green", fillOpacity = 0.7, weight = 1, group = "editable")

# mesh of possible camera locations
# different spacing (1km, 2km, and 4km
potential_cams_mesh_1km <- make.mask(type = "polygon", poly = potential_cams, spacing = 1000)
potential_cams_mesh_SigmaF <- make.mask(type = "polygon", poly = potential_cams, spacing = 2000)
potential_cams_mesh_SigmaM <- make.mask(type = "polygon", poly = potential_cams, spacing = 4000)

#buffered versions
potential_cams_mesh_SigmaF.buff <- make.mask(type = "polybuffer", poly = potential_cams, spacing = 2000, buffer = 2000)

plot(potential_cams_mesh_1km)

save(potential_cams_mesh_1km, potential_cams_mesh_SigmaF, potential_cams_mesh_SigmaM, file = "potential_cams_mesh.Rdata")

# plot on map
potential_cams_mesh_sf <- st_as_sf(potential_cams_mesh_SigmaM, coords = c("x", "y"), crs = my_crs)
leaflet()%>%
  addProviderTiles(provider= "CartoDB.Positron") %>%
  addPolygons(data = potential_cams %>% st_transform(4326), color = "red", fill = NA, weight = 2) %>%
  addCircles(data = potential_cams_mesh_sf %>% st_transform(4326), color = "green", fillOpacity = 0.7, weight = 1, group = "editable")

########################################################
#create simple CLT polygon
x <- readRDS("CLT/clt_surveyregion_anita.Rds")

# convert to UTM so can take intersections, unions etc
my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"
x1_utm <- x[1,] %>% st_transform(my_crs)
x2_utm <- x[2,] %>% st_transform(my_crs)
x3_utm <- x[3,] %>% st_transform(my_crs)
x4_utm <- x[4,] %>% st_transform(my_crs)
x5_utm <- x[5,] %>% st_transform(my_crs)

#seems to have the CRS projection but make_grid gives an error
CLTSimple <- x[1,]

#this now works, think it has CRS detail?
CLTSimple  <- x1_utm

#complex, realistic CLT polygon
#Grabouw removed but have left other small areas
#after st_transform seems to lose the CRS
CLTComplex <- x1_utm %>% 
  st_difference(x2_utm) 

# mesh of possible camera locations required for GA algorithm
# different spacing (2/3 sigma for sigmas of 2/5K and the avg) 
# spacing for grid designs used with CLT area polygon
CLTsimp_trap_locs_SigmaF <- make.mask(type = "polygon", poly = CLTSimple, spacing = 2/3 *(2000), buffer = 0)
CLTsimp_trap_locs_SigmaM <- make.mask(type = "polygon", poly = CLTSimple, spacing = 2/3 *(5000), buffer = 0)
CLTsimp_trap_locs_Avg <- make.mask(type = "polygon", poly = CLTSimple, spacing = 2/3 *(3500), buffer = 0)
CLTsimp_trap_locs_1500 <- make.mask(type = "polygon", poly = CLTSimple, spacing = 1500, buffer = 0)

#buffered versions for buffered masks
CLTsimp_trap_locs_SigmaF.buff <- make.mask(type = "polybuffer", poly = CLTSimple, spacing = 2/3 *(2000), buffer = 3*2000)
CLTsimp_trap_locs_SigmaM.buff <- make.mask(type = "polybuffer", poly = CLTSimple, spacing = 2/3 *(5000), buffer = 3*5000)
CLTsimp_trap_locs_Avg.buff <- make.mask(type = "polybuffer", poly = CLTSimple, spacing = 2/3 *(3500), buffer = 3*3500)
CLTsimp.buff <- make.mask(type = "polybuffer", poly = CLTSimple, spacing = 1000, buffer = 12000)
CLTcomp.buff <- make.mask(type = "polybuffer", poly = CLTComplex, spacing = 1000, buffer = 12000)

#one fine mesh for integration, for simple and complex
CLTsimp.nobuff <- make.mask(type = "polybuffer", poly = CLTSimple, spacing = 1000, buffer = 0)
CLTcomp.nobuff <- make.mask(type = "polybuffer", poly = CLTComplex, spacing = 1000, buffer = 0)

#convert from mask to trap objects as trap objs required by GAoptim
CLTsimp_trap_locs_SigmaF <- read.traps(data = data.frame(CLTsimp_trap_locs_SigmaF), detector = "count")
CLTsimp_trap_locs_SigmaM <- read.traps(data = data.frame(CLTsimp_trap_locs_SigmaM), detector = "count")
CLTsimp_trap_locs_Avg <- read.traps(data = data.frame(CLTsimp_trap_locs_Avg), detector = "count")

save(CLTSimple, CLTsimp.buff, CLTsimp.nobuff, CLTComplex, CLTcomp.buff, CLTcomp.nobuff, CLTsimp_trap_locs_SigmaF, CLTsimp_trap_locs_SigmaM, CLTsimp_trap_locs_Avg,
     CLTsimp_trap_locs_SigmaF.buff, CLTsimp_trap_locs_SigmaM.buff, CLTsimp_trap_locs_Avg.buff,
     file = "CLTAreas.RData")


#areas seperately
my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"

N <- read_sf(dsn = "CLT/SeperateAreas", layer = "North Survey")
N_utm <- N[1,] %>% st_transform(my_crs)

S <- read_sf(dsn = "CLT/SeperateAreas", layer = "South Survey incl unsuitable")

S <- read_sf(dsn = "CLT/SeperateAreas", layer = "South Survey minus unsuitable")
plot(st_geometry(S), axes= TRUE)

S_utm <- S[1,] %>% st_transform(my_crs)



Simp <- read_sf(dsn = "CLT/SeperateAreas", layer = "Unsuitable")
plot(st_geometry(Simp), axes= TRUE)


# mesh of possible camera locations required for GA algorithm
# spacing of 1.5 km
CLT_S_trap_locs <- make.mask(type = "polygon", poly = S_utm, spacing = 1500, buffer = 0)

CLT_N_trap_locs <- make.mask(type = "polygon", poly = N_utm, spacing = 1500, buffer = 0)

#create fine mesh for integration
CLTsimp.nobuff <- make.mask(type = "polybuffer", poly = Simp_utm, spacing = 1000, buffer = 0)

#convert from mask to trap objects as trap objs required by GAoptim
CLTsimp_trap_locs <- read.traps(data = data.frame(CLTsimp_trap_locs), detector = "count")

save(Simp_utm, CLTsimp_trap_locs, CLTsimp.nobuff, file = "CLT/SepAreas.RData")
