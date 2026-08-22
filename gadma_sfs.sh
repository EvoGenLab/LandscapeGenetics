# run from within conda env easySFS

./easySFS.py -i yvb_stairwayplot.vcf.gz -p pop_file.txt --preview -a > projvals.txt
# projvals.txt was imported and maximum was found in R script `filtering_2024-2-23.R`
# Values that maximize # segregating sites:
# Samples=188, values=21012
./easySFS.py -i yvb_stairwayplot.vcf.gz -p pop_file.txt -a --proj 188 --prefix yvb_run1

# for yvb_gone.recode.vcf, which we shall use for gadma
./easySFS.py -i yvb_gadma.recode.vcf -p popfile_yvb_gadma.txt --preview -a > projvals.txt
# Assessed in R. Preview says:
# Samples=146, segsites=20451
./easySFS.py -i yvb_gadma.recode.vcf -p popfile_yvb_gadma.txt -a --proj 146 --prefix yvb_gadma -f