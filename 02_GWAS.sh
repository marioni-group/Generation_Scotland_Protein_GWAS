
####################################################
####### STEP 01. PREPARE GRM - COMMAND LINE  #######
####################################################

gcta-1.94.1 --mbfile gs_mbfile.txt --make-grm --sparse-cutoff 0.05 --threads 10 --out mass_spec --keep samps.txt


####################################################
########## STEP 02. fastGWA Models - Loop  #########
####################################################

for i in $(ls Phenotypes/ -1 | sed -e 's/\.phen.*$//')
do
gcta-1.94.1 --mbfile gs_mbfile.txt \
--grm-sparse mass_spec \
--fastGWA-mlm \
--maf 0.01 \
--qcovar quant.txt \
--covar fact.txt \
--pheno ./Phenotypes/${i}.phen \
--threads 10 \
--out ./Outputs/${i}
done 
