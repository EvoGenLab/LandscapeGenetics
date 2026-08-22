# run from within conda env easySFS
# Downloaded from: https://github.com/isaacovercast/easySFS

# I believe yvb_gadma is just the base dataset with individuals from SIM & TL pops removed, then filtered for MAC=1
# The popfile said all individuals were in the same population


./easySFS.py -i yvb_gadma.recode.vcf -p popfile_yvb_gadma.txt --preview -a > projvals.txt
# Imported projvals.txt into R and found values that maximize # segregating sites
# Samples=146, segsites=20451
./easySFS.py -i yvb_gadma.recode.vcf -p popfile_yvb_gadma.txt -a --proj 146 --prefix yvb_gadma -f
