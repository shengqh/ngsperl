require(GenomicRanges)
require(tidyverse)
require(methylKit)
require(ggpubr)
require(Cairo)

params <- read.table(parSampleFile3, sep = "\t", header = F)
params %>% column_to_rownames("V2") %>% t() %>% data.frame()
project <- params$task_name
diff <- as.numeric(params$diff)
qvalue <- as.numeric(params$qvalue)
ncore <- as.numeric(params$ncore)

comparison <- read.table(parSampleFile1, sep = "\t", header = F)
groups <- read.table(parSampleFile2, sep = "\t", header = F)

comps <- comparison[1, "V2"]
grp1 <- comparison[1, "V1"]
grp2 <- comparison[2, "V1"]

groups <- groups[groups$V2 %in% c(grp1, grp2), ]
samples1 <- groups[groups$V2 == grp1, "V1"]
samples2 <- groups[groups$V2 == grp2, "V1"]
samples <- c(samples1, samples2)

#readin
input_path <- paste("..", "..", "..", "methylkitcorr", "result", sep="/")
input_file <- list.files(path = input_path, pattern = ".filtered.cpg.meth.rds$", all.files = FALSE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE)
cpg.all <- readRDS(file = paste0(input_path, "/", input_file))

#meta <- read.table(paste0("../../", project, "_meta.tsv"), sep = "\t", header = T)
#meta <- meta[order(row.names(meta)),]
treatment <- as.numeric(factor(groups[,"V2"]), levels = c(grp1, grp2)) - 1

sub_obj <- reorganize(cpg.all,
                      sample.ids = samples,
                      treatment = treatment)

sub_diff <- calculateDiffMeth(sub_obj,
                              overdispersion="MN",
                              adjust = "BH",
                              test = "fast.fisher",
                              mc.cores = ncore)

diff_res <- getMethylDiff(sub_diff, difference = diff, qvalue = qvalue, type = "all")
diff_res$direction <- ifelse(diff_res$meth.diff > 0, paste0("hypo_in_", grp1), paste0("hypo_in_", grp2))
saveRDS(sub_diff, file = paste0(comps, "_test.rds"))
saveRDS(diff_res, file = paste0(comps, "_methyldiff.rds"))

diff_res_grp2 <- getMethylDiff(sub_diff, difference = diff, qvalue = qvalue, type = "hypo")
diff_res_grp2$direction <- paste0("hypo_in_", grp2)
diff_res_grp1 <- getMethylDiff(sub_diff, difference = diff, qvalue = qvalue, type = "hyper")
diff_res_grp1$direction <- paste0("hypo_in_", grp1)
write.table(diff_res_grp2, file = paste0(comps, "_", grp2, ".dmcpgs"), sep = "\t", quote = F, row.names = F)
write.table(diff_res_grp1, file = paste0(comps, "_", grp1, ".dmcpgs"), sep = "\t", quote = F, row.names = F)


