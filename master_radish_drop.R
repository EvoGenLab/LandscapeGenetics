# Run radish on Hydra as a function of vegetation only 
# w/ generalized wishart regression
# Edwards distance
# written by Johan Ludot

# Modified by Heather Sibley to run as a job array
# and see if coefficients will stabilize.
# Many thanks to Matthew Kweskin for technical assistance.
# Syntax: Rscript --vanilla master.R $SGE_TASK_ID layer.tif droppop

# When not a job array: Rscript --vanilla master_radish.R layer.tif droppop

library(radish)
library(raster)
library(tidyverse)
library(adegenet)
library(poppr)

args=commandArgs(trailingOnly=TRUE)
#args=list("veg5_merge78.tif", "KB", 3, "radish_drop_KB/")

# Load genetic data

sites_csv <- read.csv("yvb_study_sites.csv")
abbs <- c("ANC", "BB", "CD", "HW", "KB", "LCK", "PU", "SIM","ST", "TL", "TMB")
abbs <- abbs[abbs != args[[2]]]
sites_csv <- sites_csv %>% filter(abbrev %in% abbs)

genind_YVB <- readRDS("genind.mac3.thin_latlong.pop3.RDS")
genind_YVB <- genind_YVB[pop=abbs]
genind_YVB <- missingno(genind_YVB, type = "loci")

# generate genetic distances
genpop_object <- genind2genpop(genind_YVB)

Dgen <- dist.genpop(genpop_object,method=2)

Dgen_mat <- as.matrix(Dgen)

# convert sites to spatial object
sites <- SpatialPoints(sites_csv[,c(4,3)],proj4string=CRS(as.character("+proj=longlat +datum=WGS84")))

# Load vegetation
vege <- raster(args[[1]])
lab <- names(vege)

aggregation_factor <- 10
vege <- aggregate(vege, fact=aggregation_factor, fun=modal)

vege_class <- cut(raster::values(vege), breaks = c(-0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 ,7.5, 8.5))
vege_cat <- 
  raster::ratify(raster::setValues(vege, as.numeric(vege_class)))

is.factor(vege_cat)

RAT <- levels(vege_cat)[[1]]
RAT$VALUE <- levels(vege) #explicitly define level names
levels(vege_cat) <- RAT

names(vege_cat) <- "layer"

# reproject sites to UTM

sites_utm48 <- spTransform(sites, crs(vege_cat))
sites_utm48 <- SpatialPoints(sites_utm48)

covs_vege <- stack(vege_cat)

# conductance surface
surface1 <- conductance_surface(covariates = covs_vege,
                                coords = sites_utm48,
                                directions = 8)


save.image(file=paste0(args[[4]], lab, "_drop_", args[[2]], "_", args[[3]], "_failsafe.RData"))

# model
model1 <- radish(Dgen_mat ~ layer,
                 data = surface1, 
                 conductance_model = radish::loglinear_conductance, 
                 measurement_model = radish::generalized_wishart,
                 nu=2335)

saveRDS(model1, paste0(args[[4]], "model_", lab, "_drop_", args[[2]], "_", args[[3]], ".RDS"))

modelnull <- radish(Dgen_mat ~ 1,
                    data = surface1, 
                    conductance_model = radish::loglinear_conductance, 
                    measurement_model = radish::generalized_wishart,
		                nu=2335)

anova(model1, modelnull) # compare models with ANOVA

fitted_conductance <- conductance(surface1, model1, quantile = 0.95)
saveRDS(fitted_conductance, paste0(args[[4]], "fitted_", lab, "_drop_", args[[2]], "_", args[[3]], ".RDS"))

saveRDS(modelnull, paste0(args[[4]], "modelnull_", lab, "_drop_", args[[2]], "_", args[[3]], ".RDS"))

#save.image(file=paste0("radish_", lab, "drop", args[[2]], ".RData"))
