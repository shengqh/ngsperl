rm(list=ls()) 
sample_name='DM_1'
outFile='DM_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20230419_scRNA_9061_mouse_byTiger/decontX/result/DM_1')

### Parameter setting end ###

source("scRNA_func.r")
library("SingleCellExperiment")
library("celda")
library("ggplot2")
library(patchwork)

sample_map = read_file_map("fileList1.txt")
raw_map = read_file_map("fileList2.txt")

read_sce<-function(countfile){
  cat("  read", countfile, "\n")
  lst = read_scrna_data(countfile)
  counts<-lst$counts
  sce <- SingleCellExperiment(assays=list(counts=counts))
  return(sce)
}

sample_name=names(sample_map)[1]
for (sample_name in names(sample_map)){
  cat(sample_name, "\n")
  countfile = sample_map[sample_name]
  raw_countfile = raw_map[sample_name]

  sce <- read_sce(countfile)
  sce.raw <- read_sce(raw_countfile)

  sce <- decontX(sce, background = sce.raw)

  df = data.frame(sce@metadata$decontX$estimates$all_cells$UMAP)
  df$contamination = unlist(sce@metadata$decontX$estimates$all_cells$contamination)

  g1<-ggplot(df, aes(DecontX_UMAP_1, DecontX_UMAP_2)) + geom_point(aes(color=contamination), size=1) + scale_color_gradient(low="lightgray", high="red") + theme_bw3()
  g2<-ggplot(df, aes(x = contamination)) + geom_histogram(aes(y = after_stat(density)), color="black", fill=NA, bins = 50) + geom_density() + theme_bw3()
  g<-g1+g2+plot_layout(ncol=2)
  png(paste0(sample_name, ".decontX.png"), width=3300, height=1500, res=300)
  print(g)
  dev.off()

  meta=colData(sce)
  saveRDS(meta, paste0(sample_name, ".decontX.meta.rds"))

  counts = round(decontXcounts(sce))
  saveRDS(counts, paste0(sample_name, ".decontX.counts.rds"))
}
