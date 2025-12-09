#Creating a polygon
x_coord <- c(0, 6000,6000,0,0)
y_coord <- c(0,0,6000,6000,0)
xym <- cbind(x_coord, y_coord)

#create a Polygon, wrap that into a Polygons object, then wrap that into a SpatialPolygons object
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)

#a SpatialPolygons object can have a CRS
proj4string(sps) = CRS("+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs")

#create sfc object needed for st_make_grid,  (or sf obj, like shape files)
sfc = st_sfc(st_polygon(list(rbind(c(0,0), c(6000,0), c(6000,6000), c(0,6000),c(0,0)))))
plot(sfc)

test = grid.design(sfc, grid.size = 100, nT = 30)
plot(sfc, axes = T)
plot(test[[1]],add=T)

#using CLT shapefile
test2 = grid.design(CLTSimple, grid.size = 2000, nT = 50)
plot(st_geometry(CLTSimple), axes = T)
plot(test2[[1]],add=T)
