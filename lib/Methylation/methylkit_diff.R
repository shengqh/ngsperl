rm(list=ls()) 
sample_name='bPlg_bPlg_vs_bPlg_unt'
outFile='bPlg_bPlg_vs_bPlg_unt'
parSampleFile1=''
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1='/nobackup/vickers_lab/projects/20251006_13927_DNAMethyl_mm10/MethylKitCorr/result/P13927.filtered.cpg.meth.rds'
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20251006_13927_DNAMethyl_mm10/MethylKitDiff/result/bPlg_bPlg_vs_bPlg_unt')

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
overdispersion <- ifelse(is.null(params$overdispersion), "MN", params$overdispersion)
test_method <- ifelse(is.null(params$test_method), "dss", params$test_method)
adjust <- ifelse(is.null(params$adjust), "BH", params$adjust) # fdr is alias of BH

cat(" Project: ", project, "\n",
    "Difference: ", difference, "\n",
    "Q-value: ", qvalue, "\n",
    "N-core: ", ncore, "\n",
    "Overdispersion: ", overdispersion, "\n",
    "Test method: ", test_method, "\n",
    "Adjust method: ", adjust, "\n")

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
cat("Reading methylation data from: ", parFile1, "\n")
cpg.all <- readRDS(file = parFile1)

cat("Reorganizing methylation data...\n")
sub_obj <- reorganize(cpg.all,
                      sample.ids = samples,
                      treatment = treatment)
rm(cpg.all)

# test_methods = c("F", "Chisq", "fast.fisher", "dss")
# for(test_method in test_methods){
if(test_method == "dss"){
  cat("Calculating differential methylation with test method: ", test_method, " and adjust method: ", adjust, "\n")
  sub_diff <- calculateDiffMethDSS(sub_obj,
                                    adjust = adjust,
                                    mc.cores = ncore)
}else{
  cat("Calculating differential methylation with test method: ", test_method, " and adjust method: ", adjust, "\n")
  sub_diff <- calculateDiffMeth(sub_obj,
                                overdispersion = overdispersion,
                                adjust = adjust,
                                test = test_method,
                                mc.cores = ncore)
}

cur_prefix <- paste0(sample_name, ".", test_method)
meth_rds = paste0(cur_prefix, ".rds")
saveRDS(sub_diff, file = meth_rds)

cat("Extracting all differential methylation results for test method", test_method, "...\n")
diff_res <- getMethylDiff(sub_diff, difference = difference, qvalue = qvalue, type = "all")
diff_res$direction <- ifelse(diff_res$meth.diff > 0, paste0("hyper_in_", treatment_group_name), paste0("hyper_in_", control_group_name))
diff_res = diff_res |> dplyr::arrange(pvalue)
write.table(diff_res, file = paste0(cur_prefix, ".dmcpgs.tsv"), sep = "\t", quote = F, row.names = F)
