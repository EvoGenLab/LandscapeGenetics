#
# YVB PCA, sPCA, and AMOVA
# With basic QC, including kinship
# And some IBD and phylogenetic networks that aren't in the paper
# 11 Sep 2023
#
# Will Brooks
####

#Modified by hcl
#2024 Feb 17
#using data filtered by Heather using SNPfilt
# Load packages ----

library(radish)
library(raster)
library(dplyr)
library(ggplot2)
library(adegenet)
library(poppr)
library(SNPRelate)
library(missMDA)
library(tidyverse)
library(interp)
library(sf)
library(vcfR)
library(ape)
library(phangorn)
library(pegas)

# Load data ----

#setwd("~/Desktop/RADseq_YVBU")
#load(file = "YVB.Rdata")
setwd("C:/Users/bigri/Desktop/RESEARCH/YVB/analysis")
load("C:/Users/bigri/Desktop/RESEARCH/YVB/analysis/YVB_dec2024_2-17-23.RData")

#bring in newly generated genind filtered in filtering_2024-2-23.R; dp5gq20, indiv max miss 0.8, snp max miss 0.5, thin = 1K, and mac = 1 or 3 ----
genind_YVB_og<-readRDS("genind.mac3.thin.RDS")
genind_YVB_og #188 indiv, 20998 loci, even with thinnign

genind_YVB_og@pop #no pop info initially
table(genind_YVB_og@pop)
str(genind_YVB_og)
rownames(genind_YVB_og@tab)

genind_YVB_test<-readRDS("genind.mac3.thin.RDS")
table(genind_YVB_test@pop) #table is empty

#match sample names to pop names

#bring in spatial elements; pop names of samples, locations of pop, and samples
pop_names <- read.csv("popdata_locinfo.csv", header=TRUE) #this one has 201 indiv, including non YVB
colnames(pop_names)[1]<-"Sample" #get rid of weird character; "Sample"   "Location" "Loc" 

pop_locs <- read.csv("yvb_study_sites.csv", header=TRUE) #19 in total, incd M'sia; has full name, latlong and abbrev

#to give pop to get sample in the genind, assoc Loc from pop_names to Sample from samples
samples<-data.frame("Sample"=rownames(genind_YVB_og@tab))


# Obtain spatial locations for populations
pops <- data.frame("pop"=as.character(unique(genind_YVB_og@pop))) 
key <- filter(pop_names, pop_names$Loc %in% pops$pop) #remove non yvb; filter pop names by what's found in genind @pop
latlon <- filter(pop_locs, pop_locs$abbrev %in% key$Loc) #from 19 loc to 15
coords <- left_join(pops, key, join_by(pop == Loc))
new_coords <- left_join(coords, latlon, join_by(pop == abbrev))
new_coords[1:5, ]

#better to merge by sample names; previous attempts caused wrong order
merged_data <- merge(x = samples, y = new_coords, by = "Sample", all.x = T)

#put into genind
genind_YVB_og@pop <- as.factor(merged_data$pop)
genind_YVB_og@other$xy <- merged_data[, c("Long", "Lat")]

table(genind_YVB_og@pop) #ST and TF combined into ST (Sungei Tengah), no more malaysian bulbuls
#ANC  BB  CD  HW  KB LCK  PE  PR  PU  SB SIM  ST  TL TMB  TT 
#11  18  22  19   3  19   1   2  21   2  16  19  17  17   1 

#write.table( pops, "pops.txt", quote = F, sep = "\t")

saveRDS(genind_YVB_og, "genind.mac3.thin_latlong.RDS")
genind_YVB_og #188, 20,998 SNPs


# combine sample and pop names and missingness info ----
# Get SNP genotype matrix with NAs preserved
geno_tab <- tab(genind_YVB_og, NA.method = "asis")

# Count missing SNPs per individual----
missing_counts <- rowSums(is.na(geno_tab))

# Proportion of missing SNPs per individual
missing_prop <- missing_counts / (2*nLoc(genind_YVB_og))

# Combine with sample + population info
YVB_og_info <- data.frame(
  Sample = indNames(genind_YVB_og),
  Population = pop(genind_YVB_og),
  Missing_SNPs = missing_counts,
  Missing_Proportion = round(missing_prop, 3)
)

write.table( YVB_og_info, "YVB_og_missingness.csv", quote = F, sep = ",")

# Estimate Ho and pi using genid and find avg per pop with at least 3 or 10 indiv ----

#Get avg Ho and plot against latitude
library(hierfstat)

#coords has  pop Longitude Latitude count

#convert ; pop size from coords
hier_YVB_og <- genind2hierfstat(genind_YVB_og)
# Compute stats
basic_stats <- basic.stats(hier_YVB_og)
Ho_locus_pop <- basic_stats$Ho
# Keep only pop with at least 3 indiv
Ho_locus_pop_filtered <- Ho_locus_pop[, coords$pop[coords$count >= 10], drop = FALSE] #drop means don't drop dimension if end up w 1 row or column
# Average Ho per population
Ho_pop_mean_filtered <- colMeans(Ho_locus_pop_filtered, na.rm = TRUE)
#make into df
Ho_pop_mean_filtered <-data.frame(pop = names(Ho_pop_mean_filtered),  Ho = Ho_pop_mean_filtered)
#then merge with coords
Ho_pop_mean_filtered <- merge(Ho_pop_mean_filtered, coords, by = "pop")

ggplot(Ho_pop_mean_filtered, aes(x = Latitude, y = Ho)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal() +
  labs(x = "Latitude", y = "Observed heterozygosity (Ho)",
       title = "Ho vs Latitude (populations ???10 indivs)")




#Do kinship using IBDking; needs vcf to GDS conversion ----
vcf_YVB_og<-read.vcfR("YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.vcf.gz")
length(colnames(vcf_YVB_og@gt)[-1]) #188; verify that this is the right file


# 1. Convert compressed VCF to GDS
snpgdsVCF2GDS(
  "YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.vcf.gz",
  "YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.gds",
  method = "biallelic.only"
)

# 2. Summarize the new GDS file
snpgdsSummary("YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.gds")

# 3. Open the GDS file
gds_YVB_og <- snpgdsOpen("YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.gds")



# run KING kinship ----
ibd_king <- snpgdsIBDKING(gds_YVB_og, autosome.only = FALSE, type = "KING-robust")
kinship_df <- snpgdsIBDSelection(ibd_king) #king-homo gives k0 and k1

#find PO and FS

# Kinship	Overall probability that a randomly chosen allele from one individual is identical by descent with an allele from the other	- 0.5: identical / duplicates
# - 0.25: 1st-degree relatives (parent-offspring, full sibs)
# - 0.125: 2nd-degree relatives
# - 0: unrelated

# IBS0	Proportion of SNPs where the pair shares 0 alleles identical by state (IBS)	- Parent-offspring: ~0
# - Parent-offspring: ~0
# - Full siblings: 0-0.25
# - Unrelated: higher values


# Parent-offspring: kinship ~0.25 & IBS0 near 0
po_pairs <- subset(kinship_df, kinship > 0.2 & kinship < 0.25 & IBS0 < 0.01)

# Full siblings: kinship ~0.25 & IBS0 > 0
fs_pairs <- subset(kinship_df, kinship > 0.2 & kinship < 0.25 & IBS0 >= 0.01 & IBS0 < 0.25)

# 4. Close the GDS file when done
snpgdsClose(gds_YVB_og)

#make sure pop labels and points match up
par(mar = c(5, 4, 4, 2))
plot(genind_YVB_og@other$xy, main = "Origingal Coordinates", xlab = "Longitude", ylab = "Latitude", pch = 20, col = 'blue', ylim = c(1.32,1.50))
text(genind_YVB_og@other$xy[,1], genind_YVB_og@other$xy[,2], labels = as.factor(genind_YVB_og@pop), pos = 3, cex = 0.8)

# Extract coordinates and population labels
coords <- as.data.frame(genind_YVB_og@other$xy)
colnames(coords) <- c("Longitude", "Latitude")
coords$pop <- as.factor(genind_YVB_og@pop)
coords <- aggregate(cbind(Longitude, Latitude) ~ pop, data = coords, FUN = mean)
coords$count<- table(genind_YVB_og@pop) 
  
#overlap pop points and texts onto map ----
country<-st_read("~/ArcGIS/Admin boundary/countries_shp/countries.shp")

pdf("pop_map.pdf", width = 8, height = 6)
ggplot(data = country) +
  geom_sf(fill = "gray", color = "gray") +
  coord_sf(xlim = c(103.6,104.1),ylim = c(1.25,1.47), expand = FALSE)+
  theme_minimal() +
  geom_point(data = coords, aes(x = Longitude, y = Latitude), color = "blue", size = 2) +  # Points
  geom_text(data = coords, aes(x = Longitude, y = Latitude, 
            label = paste0(pop, " (", count, ")")), 
            position = position_nudge(y = 0.01), size = 5) +  # Text labels
  ggtitle("")
dev.off()


#QC the merge
#combined<-cbind(genind_YVB_og@tab[,1:4], genind_YVB_og@pop, genind_YVB_og@other$xy)
#write.table( combined, "samples_pops_latlong.txt", quote = F, sep = "\t")

# PCA ----
#using mac 3, 188 indiv, incld pop with 3 or fewer indiv.  _80 dataset has 90% missingness filtering
#genind_YVB_PCA <- tab(genind_YVB_og, freq=TRUE, NA.method="mean")

#this uses genind_YVB_80 from below where SNP have been filtered to remove those with missing above 0.1; #188 5,568;
genind_YVB_PCA <- tab(genind_YVB_80, freq=TRUE, NA.method="mean") #11136 rows, which means 5568 SNPs?
class(genind_YVB_PCA) #"matrix" "array" 


#compare mine using imputePCA
genind_YVB_PCA[1:10,1:10] #using tab and freq method makes everything as prop, which is half values of my matrix
genind_YVB_PCA_imp[[1]][1:10,1:10]
#
pca.test <- dudi.pca(df = genind_YVB_PCA, center = TRUE, scale = T) # scannf = FALSE, nf = 5

pca.test$eig
(pca.test$eig[1:3]*100)/sum(pca.test$eig)
  
  
s.label(pca.test$li) #just square labels

#plot with screeplot
pdf("pca_plot.pdf", width = 7, height = 6)
par(mar = c(2, 1, 1, 4))
# I think this works for MANUSCRIPT
s.class(pca.test$li, fac=pop(genind_YVB_og), col=hue_pal(c(0, 360))(15),
        cpoint=.5,
        clabel=0.7,
        cellipse=0, xax=1,yax=2, addaxes = T,  grid = FALSE,  xlim = c(-100, 20),  ylim = c(-120, 60))
        
axis(side = 1, at = pretty(pca.test$li[,1]), labels = TRUE,  pos = -100 , outer=T)
axis(side = 2, at = pretty(pca.test$li[,2]), labels = TRUE, pos = -100)
mtext("PC1", side = 1, line=3)  # x-axis label
mtext("PC2", side = 2, line=-2)  # y-axis label

#box()
add.scatter.eig(pca.test$eig[1:20],inset = c(0.25,0.15), ratio = 0.2, nf=2,xax=1,yax=2)
dev.off()

#some outlier birds are migrants from Malaysia?
s.label(pca.test$li , boxes = F , addaxes = T )

#focus on the middle portion. Gross. I think this works? MANUSCRIPT
s.class(pca.test$li, fac=pop(genind_YVB_og), col=hue_pal(c(0, 360))(15), xlim = c(-20,20), ylim = c(-20,20),
        grid=FALSE, addaxes=FALSE,
        cpoint=.5,
        clabel=0.5,
        cellipse=0)
axis(side = 1, at = pretty(pca.test$li$li[,1]), labels = TRUE)
axis(side = 2, at = pretty(pca.test$li$li[,2]), labels = TRUE)


#add the middle portion as an inset to PCA ----
pdf("pca_plot_inset.pdf", width = 7, height = 6)
# Main plot
par(fig=c(0, 1, 0, 1), new=FALSE) #reset fig
par(mar = c(2, 1, 1, 4))  # Set the margins
# MANUSCRIPT
s.class(pca.test$li, fac=pop(genind_YVB_og), col=hue_pal(c(0, 360))(15),
        cpoint=.5,
        clabel=0.7,
        cellipse=0, xax=1, yax=2, addaxes=F, grid=FALSE, xlim=c(-100, 20), ylim=c(-120, 60))

# Add axes and labels for the main plot
axis(side=1, at=pretty(pca.test$li[,1]), labels=TRUE, pos=-100, outer=T)
axis(side=2, at=pretty(pca.test$li[,2]), labels=TRUE, pos=-100)
mtext("PC1", side=1, line=1)  # x-axis label
mtext("PC2", side=2, line=-5)  # y-axis label

# Add inset plot
par(fig=c(0.3, 0.6, 0.2, 0.5), new=TRUE)  # Adjust inset position (left, right, bottom, top)
par(mar=c(2, 1, 1, 1))  # Set margins for the inset

s.class(pca.test$li, fac=pop(genind_YVB_og), col=hue_pal(c(0, 360))(15), 
        cpoint=.5, 
        clabel=0.5, 
        cellipse=0, addaxes=F, grid=FALSE, xlim=c(-15, 15), ylim=c(-15, 15))
dev.off()

#try my method; same,? need to compare agains using the tab method to fill in missing data
#genind_YVB_PCA_imp<-imputePCA(genind_YVB_og@tab,scale=F) #this gives a list
#class(genind_YVB_PCA_imp) #list

#yvb.pca<-dudi.pca(genind_YVB_PCA_imp[[1]],scale=T,center=T, nf = 5)

#extract as axis
yvb.axes<-pca.test$li

#a cleaner plot
ggplot(yvb.axes, aes(x=Axis1 ,y=Axis2)) + geom_point(shape = 19)  +
  geom_text(aes(label=rownames(yvb.axes),size=3),color="black",size=3, vjust=-0.1)

#some outlier birds are migrants from Malaysia?
s.label(pca.test$li , boxes = F , addaxes = T )

#even more focused
s.class(pca.test$li, fac=pop(genind_YVB_og), col=funky(19), xlim = c(-25,25), ylim = c(-25,25),
        cpoint=.5,
        clabel=0.5,
        cellipse=0)
axis(side = 1, at = pretty(pca.test$li[,1]), labels = TRUE)
axis(side = 2, at = pretty(pca.test$li[,2]), labels = TRUE)


# sPCA ----
#impute missing data
#yvb_og is mac mac=3 and thin; why cutoff = 0.1 why so stringent?
genind_YVB_80 <- missingno(genind_YVB_og, cutoff = 0.1, type = "loci") #max .1 missing allowed; removes all loci containing > 0.1 missing data in the entire data set.
(genind_YVB_80) #188 5,568; why is sPCA more stringent?
genind_YVB_80@tab <- tab(genind_YVB_80, freq=TRUE, NA.method="mean")
genind_YVB_80@tab[,1:10]
genind_YVB_jitter<-genind_YVB_80

#conduct another pca with this 5568 snp data, see PCA codes above


#add jittering, spca cannot have duplicate locations
genind_YVB_jitter@other$xy[, 1] <- jitter(genind_YVB_80@other$xy[, 1], amount = 0.003)#amt of 0.003 doesn't cause duplicate loc in spca
genind_YVB_jitter@other$xy[, 2] <- jitter(genind_YVB_80@other$xy[, 2], amount = 0.003)
genind_YVB_jitter

#the jitter has 5568 loci, 188 indiv
saveRDS(genind_YVB_jitter, "genind.mac3.thin_latlong_jitter.RDS")

#visualize jitter
plot(genind_YVB_jitter@other$xy, main = "Jittered Coordinates", xlab = "Longitude", ylab = "Latitude", pch = 20, col = 'blue')

#run: adegenetTutorial(which="spca")
#how to bring lat long into genind
#Alternatively, these coordinates can be stored inside the
#genind/genpop object, preferably as @other$xy, in which case the spca function will detect
#and use this information, and not request an xy argument.
#cannot have duplicate location

#nfposi and nfnega just says how many glocal and local axes kept and retained in the object.
#mySpca <- spca(genind_YVB_jitter, ask=T, type=1, scannf=FALSE, nfposi=1, nfnega=1)
#mySpca2 <- spca(genind_YVB_jitter, ask=F, type=1, scannf=FALSE, nfposi=2, nfnega=2)

#just use mySpca3 since it retains more info
mySpca3 <- spca(genind_YVB_jitter, ask=F, type=1, scannf=FALSE, nfposi=3, nfnega=3)

#should i do just 1 global axis? Since local structure is not sig?
#mySpca4 <- spca(genind_YVB_jitter, ask=F, type=1, scannf=FALSE, nfposi=1, nfnega=0)

barplot(mySpca$eig, main="sPCA eigenvalues", col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

barplot(mySpca3$eig, main="sPCA eigenvalues", col=spectral(length(mySpca3$eig)))
legend("topright", fill=spectral(2),
       leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")

#screeplot(mySpca) #don't run - hangs

mySpca
mySpca3

head(mySpca$eig)
length(mySpca$eig) #187, 188 indiv minus 1; diff eigen values depending on if I keep 3 or 1 axes
(mySpca$c1)[1:10,] #loading of loci
(mySpca$li)[1:10,] #scores of indiv/entities
(mySpca3$li)[1:10,] #retain 3.
mySpca$as[1:10,] #pca axes onto spca axes
head(mySpca$ls) #lag vectors of the scores can be used to better perceive global structures. Lag vectors
head(mySpca$lw)

#test of global and local patterns
Gtest<-global.rtest(mySpca$tab,mySpca$lw, nperm = 999 ) #lw is list of spatial weights
Gtest
plot(Gtest, breaks = 40) #highly significant

Ltest<-local.rtest(mySpca$tab,mySpca$lw, nperm = 999 )
Ltest
plot(Ltest) #not sigficant

# plot(mySpca) don't plot, will try to create screeplot

#for myspca3; num of lw is 188 no matter what
Gtest3<-global.rtest(mySpca3$tab,mySpca3$lw, nperm = 999 ) #lw is list of spatial weights
Gtest3
plot(Gtest3, breaks = 40)

Ltest3<-local.rtest(mySpca3$tab,mySpca3$lw, nperm = 999 )
Ltest3
plot(Ltest3)

# more jitter added in plot, just a bit in analysis 
#cex =3 shows 3 scores. Global?
colorplot(mySpca, cex = 3) #only has 1 global, 1 local
colorplot(other(genind_YVB_jitter)$xy, mySpca$li, 1, cex = 1, main="colorplot of mySpca, first global score") #as above

#global, then local, retained 6, 3 global, 3 local
colorplot(other(genind_YVB_jitter)$xy, mySpca3$li, 1:3, cex = 3) #max of 3 columns
colorplot(other(genind_YVB_jitter)$xy, mySpca3$li, 4:6, cex = 3)

#colorplot(jitter(xy, factor=0.5), mySpca$ls, axes=1:2, transp=T, alpha=0.75, cex=3, add=F, main="Colorplot of PC (lagged scores)")

#interpolation; just the xy coordinates from other slot
x<-other(genind_YVB_jitter)$xy[,1]
y<-other(genind_YVB_jitter)$xy[,2]

temp <- interp(x, y, mySpca3$li[,1]) #first global axis
image(temp, col=azur(100))
#contour_levels <- c(-4, -3, -2, -1, 0, 1, 2)
contour_levels <- c( -1, 0, 1)
contour(temp, levels = contour_levels, add= T, col = 'red', labcex = 1.5)
#contour(temp, levels = contour_levels, add= F, col = 'blue')
points(x,y, pch = 16)

temp <- interp(x, y, mySpca3$li[,2]) #2nd global axis
image(temp, col=azur(100))
points(x,y)

## spca plot -----
#make very fine interpolation, but too fine
interpX <- seq(min(x),max(x),le=50) #le is length.out
interpY <- seq(min(y),max(y),le=50)
#temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY) #ls is the lagged PC, xo : sequence of x locations for rectangular output grid
temp <- interp(x, y, mySpca3$li[,1], xo=interpX, yo=interpY) #ls is the lagged PC, xo : sequence of x locations for rectangular output grid

# MANUSCRIPT
pdf("spca_interpolation.pdf", width = 6, height = 4)
par(mar = c(4.1, 4.1, 1, 2.1))
image(temp, col=colorRampPalette(c("yellow", "blue"))(100), xlab = "Longitude", ylab = "Latitude")
contour_levels <- c( -1, 0, 1)
contour(temp, levels = contour_levels, add= T, col = 'red', labcex = 1.5)
points(x,y, pch = 16, cex = 0.5)
dev.off()

#Network analysis

yvb.dist<-dist.gene(as.matrix(genind_YVB_og@tab), method = "per")
#mala.dist<-dist.gene(as.matrix(mala.gt))
#yvb.dist2<-dist(as.matrix(genind_YVB_80@tab))
yvb.nnet<-neighborNet(yvb.dist)

writeDist(yvb.dist, "yvb_dist.dist", format = "phylip")


#create a multiphylo for consensus net analysis ----
# Load required libraries
library(adegenet)
library(ape)

# Number of bootstrap replicates
nboot <- 100

# Storage for trees
trees <- vector("list", nboot)

# Number of loci
n_loci <- 5000

for (i in 1:nboot) {
  # Step 1: Resample loci (with replacement)
  sampled_loci <- sample(1:n_loci, n_loci, replace = TRUE)
  
  # Step 2: Subset genind object by loci
  gen_boot <- genind_YVB_og[, sampled_loci]
  
  # Step 3: Compute distance matrix (Euclidean distance is commonly used)
  dist_matrix <- dist.gene(as.matrix(gen_boot@tab), pairwise.deletion = T)
  
  # Step 4: Build NJ tree
  tree <- nj(dist_matrix)
  
  # Step 5: Store tree
  trees[[i]] <- tree
}

# Convert list of trees to multiPhylo object
class(trees) <- "multiPhylo"

yvb.cnet<-consensusNet(trees, prob = 0.9)

plot(yvb.cnet, show.tip.label = F, type = "equal angle", cex = 0.3)


#find out contributions of SNPs; nothing super stand out
myLoadings <- mySpca3$c1[,1]^2 #sq loadings of SNPs
names(myLoadings) <- rownames(mySpca3$c1)
loadingplot(myLoadings, xlab="Alleles",
            ylab="Weight of the alleles",
            main="Contribution of alleles \n to the first sPCA axis")

#show kernel dist; not working; needs a plot first like image
s.kde2d(genind_YVB_jitter@other$xy, add.plot = T)
points(genind_YVB_jitter@other$xy, col="red",pch=20)


#run basic IBD between indiv ----
Dgen<-dist(genind_YVB_og) #what is the distance measure? Euceldian? allele sharing
#Dgen<-dist(genind_YVB_PCA)#try PCA with interpolation

#with jittering
Dgeo<-dist(genind_YVB_jitter@other$xy) #with jittering
ibd <- mantel.randtest(Dgen,Dgeo)
ibd #Simulated p-value: 0.527 
plot(ibd)

#without jittering
Dgeo<-dist(genind_YVB_og@other$xy) #without jittering
ibd <- mantel.randtest(Dgen,Dgeo)
ibd #Simulated p-value: 0.549 
plot(ibd)

plot(Dgeo, Dgen) #within sites has much lower gen distances

#Distance between populations; filter for at least 10 indiv; also not significant ----
pop_sizes <- table(pop(genind_YVB_og))
pops_to_keep <- names(pop_sizes[pop_sizes >= 10])

genind_YVB_og_10 <- genind_YVB_og[pop(genind_YVB_og) %in% pops_to_keep, ]

genpop_YVB_og_10<-genind2genpop(genind_YVB_og_10)
Dgenpop_10<-dist.genpop(genpop_YVB_og_10)
Dgenpop_10

coords_10<-coords[coords$count >= 10, ]
Dgeopop_10<-dist(coords_10[,c("Longitude", "Latitude")])
Dgeopop_10

plot(Dgeopop_10, Dgenpop_10)
ibdpop <- mantel.randtest(Dgeopop_10, Dgenpop_10)
ibdpop #Simulated p-value: 0.32
plot(ibdpop)

# test HWE ----
#I need to find out HWE and Fist of each population, with all loci combined
#Should I use mac = 1 or mac =3?
#drop some pop?

# full population per locus

HWE.test.all <- round(pegas::hw.test(genind_YVB_og, B = 0), digits = 3)
subset(HWE.test.all,HWE.test.all[,3]<0.05)

#proportion out of HWE

nrow(subset(HWE.test.all,HWE.test.all[,3]<0.05))/nrow(HWE.test.all)

# Chi-squared test: p-value

HWE.test <- data.frame(sapply(seppop(genind_YVB_og), 
                              function(ls) pegas::hw.test(ls, B=0)[,3]))

HWE.test.chisq <- t(data.matrix(HWE.test))

#proportion of pops out of HWE per locus
alpha=0.05
Prop.loci.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 2, mean,na.rm=TRUE))
Prop.loci.out.of.HWE             

#proportion of loci out of HWE per population
Prop.pops.out.of.HWE <- data.frame(Chisq=apply(HWE.test.chisq<alpha, 1, mean,na.rm=TRUE))
Prop.pops.out.of.HWE     

# Hierarchical AMOVA ----
#not done yet
amova.result <- poppr::poppr.amova(genind_YVB, hier = ~site, 
                                   within=FALSE, method = "ade4")
amova.result

amova.test <- ade4::randtest(amova.result, nrepet = 999)
amova.test
