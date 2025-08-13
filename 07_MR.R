

# load required libraries
library(TwoSampleMR)
library(reshape)

# read in list of outcomes and their file names
outcome_list <- read.table("outcome_list.txt", header=TRUE, sep="\t")

# replace spaces in outcome names with underscores
outcome_list$Phenotype <- gsub(" ", "_", outcome_list$Phenotype)

# read in file with sample sizes for sumstats that don't have an N column
samplesizes <- read.table("samplesizes.txt", header=TRUE)

# read in exposure file
exp_dat <- read_exposure_data(
    filename = "snps_for_publication.csv",
    sep = ",",
    snp_col = "rsid",
    beta_col = "bJ",
    se_col = "bJ_se",
    effect_allele_col = "refA",
    eaf_col = "freq",
    pval_col = "pJ",
    gene_col = "rsid_Gene",
    phenotype_col = "Trait",
    chr_col = "Chr"
)

# create empty tables that will get filled with MR results
all_res <- data.frame(id.exposure=character(), id.outcome=character(), outcome=character(), exposure=character(), method=character(), nsnp=numeric(), b=numeric(), se=numeric(), pval=numeric())
all_egger <- data.frame(id.exposure=character(), id.outcome=character(), outcome=character(), exposure=character(), egger_intercept=numeric(), se=numeric(), pval=numeric())
all_steiger <- data.frame(id.exposure=character(), id.outcome=character(), outcome=character(), exposure=character(), snp_r2.exposure=numeric(), snp_r2.outcome=numeric(), correct_causal_direction=logical(), steiger_pval=numeric())

# outcome sumstats mostly in different formats to each other
# those with matching column names are grouped and run in loops and what remains are run individually

# extract group of files based on matching file name suffix
outcome_group <- outcome_list[grep("\\.f.tsv.gz", outcome_list$file),]

# start loop
for(i in 1:nrow(outcome_group)){
	# read in outcome file
	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
	)

	# add name of outcome to table
	outcome_dat$outcome <- outcome_group$Phenotype[i]
	# add sample size to table as missing in sumstats
	outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

	# harmonise outcome and exposure data
	# if there is no EAF info in sumstats then specify to use harmonising method that doesn't use EAF
	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	# run MR
	res <- mr(dat)
	all_res <- rbind(all_res, res)

	# horizontal pleiotropy test
	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	# run steiger directionality test
	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- outcome_list[grep("BUN|CKD|eGFR|urate|UACR", outcome_list$file),]

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("Disease_sum_stats/", outcome_group$file[i]),
    sep = " ",
    snp_col = "RSID",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value",
    samplesize_col = "n_total_sum"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- outcome_list[grep("combined_1000G_density_formatted_21-03-29", outcome_list$file),]

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "a2",
    other_allele_col = "a1",
    pval_col = "p-value",
    samplesize_col = "n"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- outcome_list[grep("whradjbmi|bmi.giant-ukbb", outcome_list$file),]

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "Tested_Allele",
    other_allele_col = "Other_Allele",
    eaf_col = "Freq_Tested_Allele",
    pval_col = "P",
    samplesize_col = "N"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}
	
	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- outcome_list[grep("UKB-ICBP", outcome_list$file),]

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P",
    samplesize_col = "TotalSampleSize"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}

outcome_group <- outcome_list[grep("20161107.txt.gz", outcome_list$file),]

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele2",
    other_allele_col = "Allele1",
    pval_col = "P.value"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]
	outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- read.table("outcome_group.txt", header=TRUE, sep="\t")

outcome_group$Phenotype <- gsub(" ", "_", outcome_group$Phenotype)

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("Disease_sum_stats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "ID",
    beta_col = "EFFECT",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "ALT_FREQ",
    pval_col = "P"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]
	outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- outcome_list[grep("jointGwasMc", outcome_list$file),]

for(i in 1:nrow(outcome_group)){

	outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("Disease_sum_stats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P-value",
    samplesize_col = "N"
	)

	outcome_dat$outcome <- outcome_group$Phenotype[i]

	if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
		dat <- harmonise_data(exp_dat, outcome_dat, action=1)
	} else {
		dat <- harmonise_data(exp_dat, outcome_dat)
	}

	res <- mr(dat)
	all_res <- rbind(all_res, res)

	egger <- mr_pleiotropy_test(dat)

	if(egger$outcome[1]==outcome_dat$outcome[1]){
		all_egger <- rbind(all_egger, egger)
	}

	steiger <- directionality_test(dat)

	if(steiger$outcome[1]==outcome_dat$outcome[1]){
		all_steiger <- rbind(all_steiger, steiger)
	}
}


outcome_group <- outcome_list[grep("ALLFX|ALLOA", outcome_list$file),]

outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("Disease_sum_stats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "MarkerName",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value"
    )

outcome_dat$outcome <- outcome_group$Phenotype[i]
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_group <- outcome_list[grep("PGC_", outcome_list$file),]

outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "OR",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P"
    )

outcome_dat$outcome <- outcome_group$Phenotype[i]
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]
outcome_dat$beta.outcome <- log(outcome_dat$beta.outcome)

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_group <- outcome_list[grep("ADHD|mdd", outcome_list$file),]

outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/", outcome_group$file[i]),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "OR",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    samplesize_col = "N"
    )

outcome_dat$outcome <- outcome_group$Phenotype[i]
outcome_dat$beta.outcome <- log(outcome_dat$beta.outcome)

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/Phelan_sumstats.txt",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "overall_OR",
    se_col = "overall_SE",
    effect_allele_col = "Effect",
    other_allele_col = "Baseline",
    eaf_col = "EAF",
    pval_col = "overall_pvalue"
    )

outcome_dat$outcome <- "Epithelial_ovarian_cancer"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}



outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = paste0("../disease_sumstats/30598549-GCST006979-EFO_0009270-build37.f.tsv.gz"),
    sep = "\t",
    snp_col = "snp.1",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
	)

outcome_dat$outcome <- "Bone_mineral_density"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
	dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
	dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/HbA1c_METAL_European.txt.gz",
    sep = "\t",
    snp_col = "snp",
    beta_col = "beta",
    se_col = "stderr",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf_hapmap_CEU",
    pval_col = "pvalue"
    )

outcome_dat$outcome <- "HbA1C"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
	dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
	dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/MAGIC_ln_fastingProinsulin.txt.gz",
    sep = "\t",
    snp_col = "snp",
    beta_col = "effect",
    se_col = "stderr",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pvalue"
    )

outcome_dat$outcome <- "Fasting_proinsulin"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
	dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
	dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/Fetal_BW_European_meta.NG2019.txt.gz",
    sep = " ",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ea",
    other_allele_col = "nea",
    eaf_col = "eaf",
    pval_col = "p",
    samplesize_col = "n"
    )

outcome_dat$outcome <- "Offspring_birth_weight"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/Maternal_BW_European_meta.NG2019.txt.gz",
    sep = " ",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "ea",
    other_allele_col = "nea",
    eaf_col = "eaf",
    pval_col = "p",
    samplesize_col = "n"
    )

outcome_dat$outcome <- "Own_birth_weight"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "Disease_sum_stats/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "Tested_Allele",
    other_allele_col = "Other_Allele",
    eaf_col = "Freq_Tested_Allele_in_HRS",
    pval_col = "P",
    samplesize_col = "N"
    )

outcome_dat$outcome <- "Height"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/cad.add.160614.website.txt",
    sep = "\t",
    snp_col = "markername",
    beta_col = "beta",
    se_col = "se_dgc",
    effect_allele_col = "effect_allele",
    other_allele_col = "noneffect_allele",
    eaf_col = "effect_allele_freq",
    pval_col = "p_dgc"
    )

outcome_dat$outcome <- "Coronary_Artery_Disease_additive_model"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/ShahS_31919418_HeartFailure.gz",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "N"
    )

outcome_dat$outcome <- "Heart_failure"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/META_STAGE1_GWASHR_SUMSTATS.txt",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "EFFECT_ALLELE",
    other_allele_col = "OTHER_ALLELE",
    eaf_col = "EAF",
    pval_col = "P_VALUE",
    samplesize_col = "N_TOTAL"
    )

outcome_dat$outcome <- "Heart_rate"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/Mahajan.NatGenet2018b.T2D.European.txt",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "EA",
    other_allele_col = "NEA",
    eaf_col = "EAF",
    pval_col = "Pvalue"
    )

outcome_dat$outcome <- "Type_2_Diabetes"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/HanY_prePMID_asthma_UKBB.txt.gz",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "OR",
    se_col = "SE",
    effect_allele_col = "EA",
    other_allele_col = "NEA",
    eaf_col = "EAF",
    pval_col = "Pvalue",
    samplesize_col = "N"
    )

outcome_dat$outcome <- "Asthma"
outcome_dat$beta.outcome <- log(outcome_dat$beta.outcome)

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "Disease_sum_stats/EAGLE_AD_GWAS_results_2015.txt.gz",
    sep = "\t",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "reference_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "p.value",
    samplesize_col = "European_N"
    )

outcome_dat$outcome <- "Eczema"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/Shrine_30804560_UKBiobank_FEV1_to_FVC_RATIO.txt.gz",
    sep = "\t",
    snp_col = "#SNP",
    beta_col = "beta",
    se_col = "SE_GC",
    effect_allele_col = "Coded",
    other_allele_col = "Non_coded",
    eaf_col = "Coded_freq",
    pval_col = "P_GC"
    )

outcome_dat$outcome <- "Lung_function_FEV1/FVC"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/MTAG_glaucoma_four_traits_summary_statistics.txt",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "Effect_allele",
    other_allele_col = "Non_Effect_allele",
    pval_col = "P"
    )

outcome_dat$outcome <- "Glaucoma_multi-trait_analysis"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/EvangelouE_31358974.txt",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "effect",
    se_col = "stderr",
    effect_allele_col = "allele1",
    other_allele_col = "allele2",
    eaf_col = "freq1",
    pval_col = "P",
    samplesize_col = "totalsamplesize"
    )

outcome_dat$outcome <- "Alcohol_consumption"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/AD_sumstats_Jansenetal_2019sept.txt.gz",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "EAF",
    pval_col = "P",
    samplesize_col = "Nsum"
    )

outcome_dat$outcome <- "Alzheimers_disease"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/Saxena_fullUKBB_Insomnia_summary_stats.txt",
    sep = " ",
    snp_col = "SNP",
    beta_col = "BETA_INSOMNIA",
    se_col = "SE_INSOMNIA",
    effect_allele_col = "ALELLE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P_INSOMNIA"
    )

outcome_dat$outcome <- "Insomnia_symptoms"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/migraine_ihgc2021_gws_gwama_0.txt",
    sep = "\t",
    snp_col = "rs_number",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "reference_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "p.value",
    samplesize_col = "n_samples"
    )

outcome_dat$outcome <- "Migraine"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}



outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/pgc-bip2021-all.vcf.tsv",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "PVAL",
    samplesize_col = "N"
    )

outcome_dat$outcome <- "Bipolar_disorder"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "Disease_sum_stats/Menopause_HapMap2_DayNG2015_18112015.txt.gz",
    sep = "\t",
    snp_col = "MarkerName",
    beta_col = "effect",
    se_col = "stderr",
    effect_allele_col = "allele1",
    other_allele_col = "allele2",
    pval_col = "p",
    eaf_col = "HapMap_eaf"
    )

outcome_dat$outcome <- "Age_at_Menopause"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/GenOMICC_EUR.tsv",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    pval_col = "P.value"
    )

outcome_dat$outcome <- "critical_COVID"
outcome_dat$samplesize.outcome <- 67214

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "Disease_sum_stats/UKBB.allergy.assoc.gz",
    sep = " ",
    snp_col = "SNP",
    beta_col = "OR",
    se_col = "SE",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P"
    )

outcome_dat$outcome <- "Allergic_diseases"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]
outcome_dat$beta.outcome <- log(outcome_dat$beta.outcome)

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/nallsEtAl2019_excluding23andMe_allVariants.tab",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "p",
    eaf_col = "freq",
    samplesize_col = "N"
    )

outcome_dat$outcome <- "Parkinsons_disease"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid.gz",
    sep = " ",
    snp_col = "RSID",
    beta_col = "Effect",
    se_col = "StdErr",
    effect_allele_col = "Allele1",
    other_allele_col = "Allele2",
    eaf_col = "Freq1",
    pval_col = "P-value",
    samplesize_col = "n_total_sum"
	)

outcome_dat$outcome <- "Urinary_albumin-to-creatinine_ratio"

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
    dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
    dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/25751624-GCST005536-EFO_0001359-Build37.f.tsv.gz",
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
	)

outcome_dat$outcome <- "Type_1_diabetes"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
	dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
	dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


outcome_dat <- read_outcome_data(
    snps = exp_dat$SNP,
    filename = "../disease_sumstats/29566793-GCST005647-EFO_0000253-build37.f.tsv",
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
	)

outcome_dat$outcome <- "Amyotrophic_lateral_sclerosis"
outcome_dat$samplesize.outcome <- samplesizes[which(samplesizes$Phenotype==outcome_dat$outcome[1]),"n"]

if(nrow(outcome_dat[is.na(outcome_dat$eaf.outcome),])==nrow(outcome_dat)){
	dat <- harmonise_data(exp_dat, outcome_dat, action=1)
} else {
	dat <- harmonise_data(exp_dat, outcome_dat)
}

res <- mr(dat)
all_res <- rbind(all_res, res)

egger <- mr_pleiotropy_test(dat)

if(egger$outcome[1]==outcome_dat$outcome[1]){
	all_egger <- rbind(all_egger, egger)
}

steiger <- directionality_test(dat)

if(steiger$outcome[1]==outcome_dat$outcome[1]){
	all_steiger <- rbind(all_steiger, steiger)
}


# rename columns in egger table
names(all_egger)[6:7] <- c("egger_intercept_se", "egger_intercept_pval")

# create empty results table
results <- data.frame(matrix(vector(), 0, 9, dimnames=list(c(), c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval"))))

# loop through each outcome
for(i in 1:nrow(outcome_list)){
	# extract results which are from specified outcome
    results_tmp <- all_res[which(all_res$outcome==outcome_list[i]),]
    results_tmp <- results_tmp[order(results_tmp$exposure),]
    # make a list of exposures tested with specified outcome
    exposure_list <- results_tmp[!duplicated(results_tmp$exposure),"exposure"]
    # loop through list of exposures
    for(j in 1:length(exposure_list)){
    	# extract results which are from specified exposure
        results_exp <- results_tmp[which(results_tmp$exposure==exposure_list[j]),]
        # if the specified exposure-outcome pair has a IVW method result and that result has a pvalue below the significance threshold
        # then add all method results of specified exposure-outcome pair to final results table
        if(nrow(results_exp[which(results_exp$method=="Inverse variance weighted"),])==1){
            if(results_exp$pval[results_exp$method=="Inverse variance weighted"]<0.00000159){
            results <- rbind(results, results_exp)
            }
        }
    }
}

# read in GWAS significant hits table
gwas_sig_hits <- read.csv("snps_for_publication.csv", header=TRUE)

# extract protein annotation from GWAS significant hit file
trait_genes <- gwas_sig_hits[,c(13,15)]
trait_genes <- trait_genes[!duplicated(trait_genes$Trait),]

# merge MR results and protein annotation
results <- merge(results, trait_genes, by.x="exposure", by.y="Trait", all.x=TRUE)
results <- results[order(results$outcome, results$exposure),]
results <- results[,c(2,3,1,10,4:9)]

# reshape MR results table to have all MR methods for one exposure-outcome pair in one row
results_wide <- reshape(results, idvar = c("id.exposure", "id.outcome", "exposure", "Gene.Names", "outcome"), timevar = "method", direction = "wide")

# rename column with number of SNPs
names(results_wide)[6] <- "nsnp"
# reformat column order of table
results_wide <- results_wide[,c(1:9,11:13,15:17,19:21,23:25)]
# remove any spaces and replace with underscores
names(results_wide)[7:21] <- gsub(" ", "_", names(results_wide)[7:21])
results_wide$outcome <- gsub(" ", "_", results_wide$outcome)

# merge MR results table with egger and steiger results
results_wide <- merge(results_wide, all_egger, by=c("outcome", "exposure"), all.x=TRUE)
results_wide <- merge(results_wide, all_steiger, by=c("outcome", "exposure"), all.x=TRUE)
results_wide <- results_wide[,c(1,2,5:21,24:26,29:32)]

# save out final MR results table
write.table(results_wide, "MR_results.txt", sep="\t", row.names=FALSE, quote=FALSE)
