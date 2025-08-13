
############################################################
################ STEP 1.0 COJO MODELS  #####################
############################################################

# Get files   
for i in *.ma
# Start loop
do
# Get file name 
A=$(echo $i | rev | cut -c 4- | rev)
# Get the chromosomes associated with the file (if significant)
chr_lines=`cat COJO_Inputs/${A}_cojo_sig_chrs.txt`
# Loop through chromosome(s) for this protein
for chr in $chr_lines
# Start sub-loop (nested loop)
do
# COJO code
gcta-1.94.1 --bfile GS20K_chr${chr}_HRC.r1-1_nomono_I4_cpra  \
--maf 0.01 \
--chr $chr \
--cojo-file COJO_Inputs/${A}.ma \
--cojo-slct \
--keep samps.txt \
--cojo-p 1.139052e-10  \
--out ../COJO_Outputs/${A}_chr${chr}
# Finish nested loop
done
# Finish overall loop
done