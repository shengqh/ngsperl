rm(list=ls()) 
sample_name='CLTI_vs_Control'
outFile='CLTI_vs_Control'
parSampleFile1=''
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1='/nobackup/brown_lab/projects/20231214_10473_Methylation_hg38/MethylKitCorr/result/P10473.filtered.cpg.meth.rds'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20231214_10473_Methylation_hg38/MethylKitDiff/result/CLTI_vs_Control')

### Parameter setting end ###

require(GenomicRanges)
require(tidyverse)
require(methylKit)
require(ggpubr)
require(Cairo)

#read parameter
params <- read.table(parSampleFile2, sep = "\t", header = F)
params <- params %>% column_to_rownames("V2") %>% t() %>% data.frame()
project <- params$task_name
difference <- as.numeric(params$difference)
qvalue <- as.numeric(params$qvalue)
ncore <- as.numeric(params$ncore)

#read comparison
comparisons <- read.table(parSampleFile3, sep = "\t", header = F)
comparison <- comparisons[comparisons$V2 == sample_name, ]
control_group_name <- comparison[1, "V1"]
treatment_group_name <- comparison[2, "V1"]

#read group
groups <- read.table(parSampleFile4, sep = "\t", header = F)
stopifnot(all(c(control_group_name, treatment_group_name) %in% groups$V2))

groups <- groups[groups$V2 %in% c(control_group_name, treatment_group_name), ]
control_names <- groups[groups$V2 == control_group_name, "V1"]
treatment_names <- groups[groups$V2 == treatment_group_name, "V1"]
samples <- c(control_names, treatment_names)
treatment <- rep(c(0, 1), c(length(control_names), length(treatment_names)))

#read meth files
cpg.all <- readRDS(file = parFile1)

sub_obj <- reorganize(cpg.all,
                      sample.ids = samples,
                      treatment = treatment)
rm(cpg.all)

sub_diff <- calculateDiffMeth(sub_obj,
                              overdispersion="MN",
                              adjust = "BH",
                              test = "fast.fisher",
                              mc.cores = ncore)
test_rds=paste0(sample_name, "_test.rds")
saveRDS(sub_diff, file = test_rds)
rm(sub_obj)
if(0){
  sub_diff <- readRDS(test_rds)
}

diff_res <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "all")
diff_res$direction <- ifelse(diff_res$meth.diff > 0, paste0("hypo_in_", control_group_name), paste0("hypo_in_", treatment_group_name))
methyldiff_rds=paste0(sample_name, "_methyldiff.rds")
saveRDS(diff_res, file = methyldiff_rds)
if(0){
  diff_res <- readRDS(methyldiff_rds)
}

diff_res_treatment_high <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "hypo")
if(nrow(diff_res_treatment_high) > 0){
  diff_res_treatment_high$direction <- paste0("hypo_in_", treatment_group_name)
  diff_res_treatment_high <- diff_res_treatment_high[order(diff_res_treatment_high$meth.diff),]
}
write.table(diff_res_treatment_high, file = paste0(sample_name, "_", treatment_group_name, ".dmcpgs"), sep = "\t", quote = F, row.names = F)

diff_res_control_high <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "hyper")
if(nrow(diff_res_control_high) > 0){
  diff_res_control_high$direction <- paste0("hypo_in_", control_group_name)
  diff_res_control_high <- diff_res_control_high[order(diff_res_control_high$meth.diff),]
}
write.table(diff_res_control_high, file = paste0(sample_name, "_", control_group_name, ".dmcpgs"), sep = "\t", quote = F, row.names = F)
