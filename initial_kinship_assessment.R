#
# YVB relatedness
# 16 May, 2024
#
# Will Brooks
############################################

setwd("~/Desktop/RADseq_YVBU/SNPrelate")

library(SNPRelate)

# find relatives with SNPrelate ----

snpgdsVCF2GDS("YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac1.thin.vcf.gz", "YVB.gds", method="biallelic.only")
#snpgdsVCF2GDS("YVB.dp5gq20.sample0.8.biallelic.snp0.5.mac3.thin.vcf.gz", "YVB.gds", method="biallelic.only")
snpgdsSummary("YVB.gds")

genofile <- snpgdsOpen("YVB.gds")

# LD prune
set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))

# Estimate IBD coefficients

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

set.seed(100)
#snp.id <- sample(snpset.id, 17789)  # subset snps if needed
ibd <- snpgdsIBDMLE(genofile, 
                    sample.id=sample.id, 
                   # snp.id=snp.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)

# Make a data.frame
ibd.coeff <- snpgdsIBDSelection(ibd)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MLE)")
lines(c(0,1), c(1,0), col="red", lty=2)

# Look at pairings with relatedness > 0.2
subset(ibd.coeff,kinship > 0.2)

# Highlight the only related pair in different sites
points(ibd.coeff[7974,]$k0, 
       ibd.coeff[7974,]$k1,
       pch=25,
       col = "black",
       bg = "red",
       cex = 1.25)

