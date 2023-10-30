rm(list=ls()) 
sample_name='B_vs_A'
outFile='B_vs_A'
parSampleFile1=''
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1='/nobackup/h_cqs/shengq2/temp/20231027_10473_WGBS_real/MethylKitCorr/result/P10473.filtered.cpg.meth.rds'
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20231027_10473_WGBS_real/MethylKitDiff/result/B_vs_A')

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
grp1 <- comparison[1, "V1"]
grp2 <- comparison[2, "V1"]

#read group
groups <- read.table(parSampleFile4, sep = "\t", header = F)
stopifnot(all(c(grp1, grp2) %in% groups$V2))

groups <- groups[groups$V2 %in% c(grp1, grp2), ]
samples1 <- groups[groups$V2 == grp1, "V1"]
samples2 <- groups[groups$V2 == grp2, "V1"]
samples <- c(samples1, samples2)

#read meth files
cpg.all <- readRDS(file = parFile1)

treatment <- as.numeric(factor(groups[,"V2"]), levels = c(grp1, grp2)) - 1

sub_obj <- reorganize(cpg.all,
                      sample.ids = samples,
                      treatment = treatment)

sub_diff <- calculateDiffMeth(sub_obj,
                              overdispersion="MN",
                              adjust = "BH",
                              test = "fast.fisher",
                              mc.cores = ncore)
saveRDS(sub_diff, file = paste0(comps, "_test.rds"))

diff_res <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "all")
diff_res$direction <- ifelse(diff_res$meth.diff > 0, paste0("hypo_in_", grp1), paste0("hypo_in_", grp2))
saveRDS(diff_res, file = paste0(comps, "_methyldiff.rds"))

diff_res_grp2 <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "hypo")
diff_res_grp2$direction <- paste0("hypo_in_", grp2)
diff_res_grp2 <- diff_res_grp2[order(diff_res_grp2$meth.diff),]
write.table(diff_res_grp2, file = paste0(comps, "_", grp2, ".dmcpgs"), sep = "\t", quote = F, row.names = F)

diff_res_grp1 <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "hyper")
diff_res_grp1$direction <- paste0("hypo_in_", grp1)
diff_res_grp1 <- diff_res_grp1[order(diff_res_grp1$meth.diff),]
write.table(diff_res_grp1, file = paste0(comps, "_", grp1, ".dmcpgs"), sep = "\t", quote = F, row.names = F)
