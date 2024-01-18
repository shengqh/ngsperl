rm(list=ls()) 
outFile='P10473'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20231109_10473_WGBS_real/MethylKitCorr/result')

### Parameter setting end ###

require(GenomicRanges)
require(tidyverse)
require(methylKit)
require(ggpubr)
require(Cairo)

params <- read.table(parSampleFile2, sep = "\t", header = F)
params <- params %>% column_to_rownames("V2") %>% t() %>% data.frame()
project <- params$task_name
assembly = params$org
mincov = as.numeric(params$mincov)

#prepare to find all specific format files
cpg.all <-read.table(parSampleFile1, sep="\t", header=F)
cpg.infile = cpg.all$V1
file.id = cpg.all$V2

meta <- read.table(parSampleFile3, sep = "\t", header = F)
colnames(meta) <- c("sample", "group")
meta <- meta[order(row.names(meta)),]
#treatment <- as.numeric(factor(meta[,var])) - 1
treatment <- rep(0, nrow(meta))
#var = params$var
var = "group"

cpg.obj <- methRead(location = as.list(cpg.infile),
                    sample.id = as.list(file.id),
                    assembly = assembly,
                    treatment = treatment,
                    context = "CpG",
                    mincov = mincov)

filtered.obj <- filterByCoverage(cpg.obj,
                                 lo.count=mincov,
                                 lo.perc=NULL,
                                 hi.count=NULL,
                                 hi.perc=99.99)
#Merging samples
filtered.cpg.meth <- methylKit::unite(filtered.obj, destrand=FALSE)
saveRDS(filtered.cpg.meth, paste0(project, ".filtered.cpg.meth.rds"))

cpg_bvalue_df <- data.frame(filtered.cpg.meth) %>%
  mutate(id = paste(chr, start, end, sep = "_"))
for(i in 1:length(file.id)){
  id <- file.id[i]
  ncs <- paste0("numCs", i)
  cov <- paste0("coverage", i)
  cpg_bvalue_df[, id] <- cpg_bvalue_df[, ncs] / cpg_bvalue_df[, cov]
}
cpg_bvalue_df <- cpg_bvalue_df %>%
  dplyr::select(one_of("id", file.id)) %>%
  column_to_rownames("id")

cpg_bvalue_cor <- cor(cpg_bvalue_df, method = "pearson")

cpg_bvalue_mds <- (1 - cpg_bvalue_cor) %>%
  cmdscale() %>%
  as_tibble()

colnames(cpg_bvalue_mds) <- c("Dim.1", "Dim.2")
#rownames(cpg_bvalue_mds) <- colnames(cpg_bvalue_df)
cpg_bvalue_mds[,var] <- meta[,var]

# Plot MDS
label=ifelse(ncol(cpg_bvalue_df) > 20, NULL, colnames(cpg_bvalue_df))
dms_plot <- ggscatter(cpg_bvalue_mds, x = "Dim.1", y = "Dim.2",
                      label = label,
                      xlab = "MDS 1",
                      ylab = "MDS 2",
                      color = var,
                      palette = "jco",
                      font.label = 8,
                      size = 1,
                      ellipse = F,
                      ellipse.type = "norm",
                      repel = TRUE) +
  theme_bw()


ggsave(paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.pdf"), dms_plot, width=4, height=3, units="in", dpi=300, bg="white")

ggsave(paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.png"), dms_plot, width=4, height=3, units="in", dpi=300, bg="white")

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
