# Run ResistanceGA on Hydra as a function of vegetation only with merged water classes
# Edwards genetic distance

# libraries

library(adegenet)
library(gdistance)
library(ResistanceGA)

args <- commandArgs(trailingOnly=TRUE)

# Load genetic data

yvb <- readRDS("genind.mac3.thin_latlong.pop3.RDS")
yvb_pop <- genind2genpop(yvb)

edwards_dist <- dist.genpop(yvb_pop, method=2)

# Load landscape data, downsample, and project to UTM

vege <- raster(args[[2]])
lab <- names(vege)
aggregation_factor <- 10
vege <- aggregate(vege, fact=aggregation_factor, fun=modal)

coords <- data.frame("pop"=yvb@pop, "Long"=yvb@other$xy$Long, "Lat"=yvb@other$xy$Lat)
coords <- unique(coords)
order <- data.frame("pop"=names(edwards_dist))
coords <- dplyr::left_join(order, coords)
coords <- SpatialPoints(coords=coords[,c("Long", "Lat")], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
crs <- "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs"
coords <- spTransform(coords, CRSobj = CRS(crs))
coords <- SpatialPoints(coords, proj4string=CRS(crs))

out <- paste("rga", lab, args[[1]], sep="_")

# set up inputs

GA.inputs <- GA.prep(vege, 
                     Results.dir=paste0(out, "/"),
                     parallel=6)
gdist.inputs <- gdist.prep(length(coords),
                           samples = coords,
                           response = as.vector(edwards_dist),
                           method="commuteDistance")

# run optimization

SS_RESULTS.gdist <- SS_optim(gdist.inputs = gdist.inputs,
                             GA.inputs = GA.inputs)
saveRDS(SS_RESULTS.gdist, paste0(out, "_results.RDS"))
#save.image(file = paste0(out, ".RData"))

