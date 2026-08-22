#
# YVB Ne
# 10 Jun, 2024
#
# Will Brooks
############################################

setwd("~/Desktop/RADseq_YVBU/NeEstimator")

library(RLDNe)
library(dplyr)

genind <- readRDS("genind.mac3.thin_latlong.RDS")
genotype <- as.data.frame(genind@tab)

# Function to change format of genotypes
transform_values <- function(x) {
  ifelse(x == 0, 11,
         ifelse(x == 1, 12,
                ifelse(x == 2, 22, x)))
}

# Apply the transformation to the entire data frame
genotype_trans <- genotype %>% mutate_all(transform_values)

# Run the NeEstimator
gp_file <- write_genepop_zlr(loci = genotype_trans,
                           pops = genind@pop,
                           ind.ids = row.names(as.data.frame(genind@tab)),
                           folder = "",
                           filename ="genepop_output.txt",
                           missingVal = NA,
                           ncode = 2,
                           diploid = T)

param_files <- NeV2_LDNe_create(input_file = gp_file$Output_File,
                               param_file = "Ne_params.txt",
                               NE_out_file = "Ne_out.txt")

run_LDNe(LDNe_params = param_files$param_file)

Ne_estimates <- readLDNe_tab(path = param_files$Ne_out_tab)

# Run again with one population

setwd("~/Desktop/RADseq_YVBU/NeEstimator_onepop")

# Run the NeEstimator
gp_file <- write_genepop_zlr(loci = genotype_trans,
                             pops = rep("singapore",length(genind@pop)),
                             ind.ids = row.names(as.data.frame(genind@tab)),
                             folder = "",
                             filename ="genepop_output.txt",
                             missingVal = NA,
                             ncode = 2,
                             diploid = T)

param_files <- NeV2_LDNe_create(input_file = gp_file$Output_File,
                                param_file = "Ne_params.txt",
                                NE_out_file = "Ne_out.txt")

run_LDNe(LDNe_params = param_files$param_file)

Ne_estimates <- readLDNe_tab(path = param_files$Ne_out_tab)

