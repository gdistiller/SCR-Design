#26 Sept 2024, 
#makes potential trap locs and masks for CLT area after getting updated mask
#does it for the two areas seperately

library(sf)
library(leaflet)
library(secr)

#areas seperately
my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"

N <- read_sf(dsn = "CLT/SeperateAreas", layer = "North Survey_edit AW")
N_utm <- N[1,] %>% st_transform(my_crs)

S <- read_sf(dsn = "CLT/SeperateAreas", layer = "South Survey minus unsuitable_edit AW")
S_utm <- S[1,] %>% st_transform(my_crs)

Both <- read_sf(dsn = "CLT/SeperateAreas", layer = "North and South_edit AW_union")
Both_utm <- Both[1,] %>% st_transform(my_crs)

#mesh of possible camera locations required for GA algorithm
#same as the mask due to a buffer of 0 
# spacing of 1.5 km needed for traps object at coarser resolution

temp  <- make.mask(type = "polygon", poly = S_utm, spacing = 1500, buffer = 0)
CLT_S_trap_locs <- read.traps(data = data.frame(temp), detector = "count")

temp  <- make.mask(type = "polygon", poly = N_utm, spacing = 1500, buffer = 0)
CLT_N_trap_locs <- read.traps(data = data.frame(temp), detector = "count")

temp  <- make.mask(type = "polygon", poly = Both_utm, spacing = 1500, buffer = 0)
CLT_Both_trap_locs <- read.traps(data = data.frame(temp), detector = "count")

CLT_S_msk <- make.mask(type = "polygon", poly = S_utm, spacing = 1000, buffer = 0)
CLT_N_msk <- make.mask(type = "polygon", poly = N_utm, spacing = 1000, buffer = 0)
CLT_Both_msk <- make.mask(type = "polygon", poly = Both_utm, spacing = 1000, buffer = 0)

CLTmsks <- list("N Polygon SF obj" = N_utm, "S Polygon SF obj" = S_utm, "Both Polygon SF obj" = Both_utm, 
                "Mask of South" = CLT_S_msk, "Mask of North" = CLT_N_msk, "Mask of Both" = CLT_Both_msk, "Trap locs for South" = CLT_S_trap_locs, "Trap locs for North" = CLT_S_trap_locs, "Trap locs for both" = CLT_Both_trap_locs)

save(CLTmsks, file = "CLT/SepAreas.RData")
