
########################################################
########## STEP 01. fastGWA-GEI Models - Loop  #########
########################################################

 for i in $(ls Phenotypes/ -1 | sed -e 's/\.phen.*$//')
do
gcta-1.94.1 --mbfile gs_mbfile.txt \
--grm-sparse mass_spec \
--fastGWA-mlm \
--noSandwich \
--maf 0.01 \
--envir envir.txt \
--qcovar quant_gei.txt \
--pheno ./Phenotypes/${i}.phen \
--threads 10 \
--out ./Sex_Interaction_Outputs/${i}
done 
