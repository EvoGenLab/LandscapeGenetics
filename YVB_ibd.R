# DEFUNCT - I believe everything in this script is covered in Johan's IBD and HC's population structure analysis and we didn't use this



# Preliminary exploration of heterozygosity and genetic structure in YVB RADseq data
# 2/2/24
#
# Heather Willis


# Load packages and setwd

library(adegenet)
library(hierfstat)
library(EcoGenetics)
library(tidyverse)
setwd("C:/Users/hrwil/OneDrive - George Mason University - O365 Production/Studies/yellow_vented_bulbul")

# Load data

bird.ind <- import2genind("pca.r50.YVB.snps.gen")
pop_names <- read.csv("popdata_locinfo.csv", header=TRUE)
pop_locs <- read.csv("yvb_study_sites.csv", header=TRUE)

# Obtain spatial locations for populations

pops <- data.frame("pop"=as.character(unique(bird.ind@pop)))
key <- filter(pop_names, pop_names$Sample %in% pops$pop)
latlon <- filter(pop_locs, pop_locs$abbrev %in% key$Loc)
coords <- left_join(pops, key, join_by(pop == Sample))
new_coords <- left_join(coords, latlon, join_by(Loc == abbrev))
new_coords


dists <- as.vector(dist(new_coords[, c("Lat", "Long")]))
dists


# Heterozygosity by population

Hobs.vec <- as.vector(readRDS("Hobs.by.pop.rds"))
mean(as.numeric(Hobs.vec))
range(Hobs.vec)
sum <- readRDS("bird.ind.summary.rds")
h <- t(data.frame(Hobs.vec))
nPop <- data.frame(sum$n.by.pop)

plot.hetzyg <- data.frame("Lat"=new_coords$Lat, "Long"=new_coords$Long, 
                          "Het"=h, "nPop"=nPop)

ggplot(data=plot.hetzyg, aes(x=Long, y=Lat, col = Het)) + 
  geom_point(size=plot.hetzyg$sum.n.by.pop) +
  theme_minimal()


## Genetic distances.

# Make the data into a genpop object
bird.pop <- genind2genpop(bird.ind)

# Genetic distances with adegenet
dist.pop.nei <- dist.genpop(bird.pop, method=1)
dist.pop.edwards <- dist.genpop(bird.pop, method=2)
dist.pop.rogers <- dist.genpop(bird.pop, method=4)


# Pairwise Fst. 
#dist.pop.pairwisefst.hierfstat <- 
# as.dist(pairwise.neifst(genind2hierfstat(bird.ind)))
#saveRDS(dist.pop.pairwisefst.hierfstat, file="dist.pop.pairwisefst.hierfstat.rds")

pairwisefst.hierfstat <- readRDS("dist.pop.pairwisefst.hierfstat.rds")

# Individual-based DPS

dps <- propShared(bird.ind)
heatmap(dps)

bird.hierfstat <- hierfstat::genind2hierfstat(bird.ind)
stats <- hierfstat::basic.stats(bird.hierfstat)
fis <- apply(stats$Fis, 2, function(x) mean(x, na.rm=TRUE))
fis

# plot pairwise Fst against physical distance
# Both outliers involve L0640, a population of 2 individuals.

fst <- as.vector(pairwisefst.hierfstat)
fst.plot <- ggplot(data.frame(dists, fst), aes(x=dists, y=fst)) +
  geom_point() + geom_smooth(col = "red") + theme_minimal() +
  xlab("Geographic Distance") + ylab("Pairwise Fst between populations")
fst.plot


# plot Nei's genetic distance against physical distance

nei <- as.vector(dist.pop.nei)
nei.plot <- ggplot(data.frame(dists, nei), aes(x=dists, y=nei)) +
  geom_point() + geom_smooth(col = "blue") + theme_minimal() +
  xlab("Geographic Distance") + ylab("Nei's genetic distance between populations")
nei.plot

# Edwards' distance
ed <- as.vector(dist.pop.edwards)
ed.plot <- ggplot(data.frame(dists, ed), aes(x=dists, y=ed)) +
  geom_point() + geom_smooth(col = "blue") + theme_minimal() +
  xlab("Geographic Distance") + ylab("Edwards' genetic distance between populations")
ed.plot


# Rogers' distance
steve <- as.vector(dist.pop.rogers)
steve.plot <- ggplot(data.frame(dists, steve), aes(x=dists, y=steve)) +
  geom_point() + geom_smooth(col = "blue") + theme_minimal() +
  xlab("Geographic Distance") + ylab("Rogers' genetic distance between populations")
steve.plot



