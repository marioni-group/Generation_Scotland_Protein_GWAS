
############################################################
######### STEP 1.0 PREPARE INPUTS - SUMMARY STATS ##########
############################################################ 

# Load requisite libraries
library(data.table)

# Extract results files 
files=list.files(".",".fastGWA")

# Define Bonferroni-threshold
threshold=(5e-8/439)

# Set up list to store all significant SNPs (pre-COJO)
list1=list()

# Loop through files and prepare them for COJO analysis 
for(i in files){ 
# Read in file   
tmp=as.data.frame(fread(i))
# Only proceed if there are significant results
tmp.test=tmp[which(tmp$P <= threshold),]
# Condition 
if(nrow(tmp.test)==0){ print("SKIPPED") } else {
# First lets store all significant results for a separate supplementary table
list1[[i]]=tmp.test
# Extract columns of interest 
tmp1 <- tmp[,c(2,4,5,7,8,9,10,6)]
# Rename columns for COJO 
names(tmp1)[4:7] <- c("freq", "b", "se", "p")
# Rename file for output 
A = gsub(".fastGWA", "", i)
# Write out file 
write.table(tmp1, paste0("../COJO_Inputs/",A,".ma"), row.names = F, quote = F)
# Print to denote completion
print(i)  
}
} 

# Combine all significant SNPs (pre-COJO)
list1=as.data.frame(do.call("rbind",list1))


################################################################
###### STEP 2.0 EXTRACT CHROMOSOMES WITH SIGNIF. SNPs ##########
################################################################

# Loop through files and find significant SNPs 
for(i in files){ 
# Read in file   
tmp=as.data.frame(fread(i))
# Extract significant associations 
tmp1 = tmp[which(tmp$P<=threshold),]
# Condition - no need to proceed if no significant results 
if(nrow(tmp1)==0){ print("SKIPPED") } else {
# Get their chromosome numbers 
chr = as.data.frame(unique(tmp1$CHR)) 
# Rename file for output 
A = gsub(".fastGWA", "", i)
# Write out file 
write.table(chr, paste0("../COJO_Inputs/", A, "_cojo_sig_chrs.txt"), row.names = F, col.names = F, quote = F)
# Print to denote completion
print(i)  
}
}

