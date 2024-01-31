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
saveRDS(sub_diff, file = paste0(sample_name, "_test.rds"))
rm(sub_obj)

diff_res <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "all")
diff_res$direction <- ifelse(diff_res$meth.diff > 0, paste0("hypo_in_", grp1), paste0("hypo_in_", grp2))
saveRDS(diff_res, file = paste0(sample_name, "_methyldiff.rds"))

diff_res_grp2 <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "hypo")
diff_res_grp2$direction <- paste0("hypo_in_", grp2)
diff_res_grp2 <- diff_res_grp2[order(diff_res_grp2$meth.diff),]
write.table(diff_res_grp2, file = paste0(sample_name, "_", grp2, ".dmcpgs"), sep = "\t", quote = F, row.names = F)

diff_res_grp1 <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "hyper")
diff_res_grp1$direction <- paste0("hypo_in_", grp1)
diff_res_grp1 <- diff_res_grp1[order(diff_res_grp1$meth.diff),]
write.table(diff_res_grp1, file = paste0(sample_name, "_", grp1, ".dmcpgs"), sep = "\t", quote = F, row.names = F)
