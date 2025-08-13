
# Load requisite libraries 
library(data.table)

# Read in phenotypes
phenos=readRDS("GS_ProteinGroups_RankTransformed_23Aug2023.rds")

####################################################################
######### STEP 01. PREPARE COVARIATE AND ENVIRONMENT FILES  ########
####################################################################

# Read in phenotypes - we want age and genetic PCs as quantitative covariates 
# Age 
cov=readRDS("GS_phenos_internal_23Aug2023_REM.rds")
cov1=cov[,c("id","age")]
# PCs
pcs=as.data.frame(fread("GS20K_ALL_MAF5_PCA.eigenvec"))
names(pcs)[1:2]=c("FID","id")
pcs$FID=NULL

# Merge information together 
cov2=merge(cov1,pcs,by="id",all.x=T)

# Subset to people with protein data 
samps=read.table("samps.txt",header=F)
cov2=cov2[which(cov2$id %in% samps$V1),]

# Create quantitative covariate file 
tmp = data.frame(fid=cov2$id, iid=cov2$id, age=cov2$age, PC1=cov2$V3, PC2=cov2$V4, PC3=cov2$V5, PC4=cov2$V6, PC5=cov2$V7, PC6=cov2$V8,
                 PC7=cov2$V9, PC8=cov2$V10, PC9=cov2$V11, PC10=cov2$V12, PC11=cov2$V13, PC12=cov2$V14, PC13=cov2$V15, PC14=cov2$V16,
                 PC15=cov2$V17, PC16=cov2$V18, PC17=cov2$V19, PC18=cov2$V20, PC19=cov2$V21, PC20=cov2$V22)


# Create gei environment file - will be used in GEI tests 
tmp1=data.frame(fid=cov$id, iid=cov$id, sex=cov$sex)
# Need to format binary variables as 0 and 1 
tmp1$sex=ifelse(tmp1$sex=="M",0,1)

# Write out files 
write.table(tmp, file="quant_gei.txt", sep='\t', row.names=F, quote=F, col.names=F)
write.table(tmp1, file="envir.txt", sep='\t', row.names=F, quote=F, col.names=F)




