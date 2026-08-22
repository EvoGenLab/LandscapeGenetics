# Filtering
# Heather Sibley 2/13/24

#use this to work on Will's data - 2/17/2024
setwd("C:/Users/bigri/Desktop/RESEARCH/YVB/analysis")
#load("C:/Users/bigri/Desktop/RESEARCH/YVB/analysis/hw_filtering.RData") #old
load("C:/Users/bigri/Desktop/RESEARCH/YVB/analysis/filtering_2024-2-23.RData") #new

library(SNPfiltR)
library(vcfR)
library(adegenet)
setwd("C:/Users/hrwil/OneDrive - George Mason University - O365 Production/Studies/yellow_vented_bulbul")


#2-23-24 import all SNPs data exported by Will
all<-read.vcfR("populations.snps.vcf")
all
#192 samples
#32901 CHROMs
#143,217 variants

# Exploratory visualizations
num_plots <- 9  # Update this with the actual number of plots
num_rows <- ceiling(num_plots / 2)
pdf("allsnps_QC.pdf", paper = "letter")

par(mfrow = c(num_rows, 2))

hard_filter(all) #some super high dp
max_depth(all)
missing_by_sample(vcfR=all)
missing_by_snp(all)

dev.off()

# Remove LCK1
names <- colnames(all@gt)
sum(c(names == "YVB_LCK1"))
names <- c(names != "YVB_LCK1")
all <- all[,names]
all
#191 samples
#32901 CHROMs
#143,217 variants

# Filter the dataset
vcf <- hard_filter(vcfR=all, depth=5, gq=20) # lower bound for read depth & genotype quality
#4.54% of genotypes fall below a read depth of 5 and were converted to NA
#0.96% of genotypes fall below a genotype quality of 20 and were converted to NA
vcf <- max_depth(vcf, maxdepth=100) # higher bound for read depth, mean cov acroos all is 34, mean+3SD (19) = 91.7 
#4.91% of SNPs were above a mean depth of 100 and were removed from the vcf
vcf 
#191 samples
#31391 CHROMs
#136,181 variants

#here are some filter combo: sample cuttoff (0.5, 0.8); snp cutoff (0.5, 0.8). mac (1,3) 
vcf <- missing_by_sample(vcfR=vcf, cutoff=0.8) # 3 samples are above a 0.8 missing data cutoff, and were removed from VCF
# Individuals filtered: PR3, TL4-7, TMB10 using 80%
# TL4-7 and TMB10 are both in large populations (18 individuals each).
# PR3 is in a population of 3. Unfortunately, PR3 has 96% missing data.

vcf <- filter_biallelic(vcf) #none removed
vcf <- missing_by_snp(vcf, cutoff=0.5) # remove snp with more than 50% missing data, cut-off, max missing data allowed
#29.53% of SNPs fell below a completeness cutoff of 0.5 and were removed from the VCF
vcf
#188 samples
#22114 CHROMs
#95,965 variants

vcf.mac1 <- min_mac(vcf, min.mac=1) # drop loci that are now monoallelic, 
vcf.mac1 #2.9% of SNPs fell below a minor allele count of 1

vcf.mac3 <- min_mac(vcf, min.mac=3) # MAC 3
vcf.mac3 #16.05% of SNPs fell below a minor allele count of 3 and were removed from the VCF

vcf.mac1.thin<-distance_thin(vcf.mac1, min.distance = 1000)
vcf.mac3.thin<-distance_thin(vcf.mac3, min.distance = 1000)
vcf.mac1.thin
#22114 out of 95965
vcf.mac3.thin
#20998 out of 80567

write.vcf(vcf.mac1.thin, file="YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac1.thin.vcf.gz")
write.vcf(vcf.mac3.thin, file="YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.vcf.gz")

genind.mac1.thin<-vcfR2genind(vcf.mac1.thin)
genind.mac3.thin<-vcfR2genind(vcf.mac3.thin)

saveRDS(genind.mac1.thin, "genind.mac1.thin.RDS")
saveRDS(genind.mac3.thin, "genind.mac3.thin.RDS")

genind.mac1.thin@pop


#see max num of NA per column
prop_na_per_col<-apply(bird.ind50@tab, 2, function(x) max(sum(is.na(x)/length(x))))
max(prop_na_per_col) #0.9947644
boxplot(prop_na_per_col, main = "Prop of NA per SNP")

names(prop_na_per_col)[which.max(prop_na_per_col)]
bird.ind50@tab[, which(colnames(bird.ind50@tab) == "69_152.01")] #has a ton of NAs

