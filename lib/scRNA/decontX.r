rm(list=ls()) 
sample_name='iSGS_3364_1'
outFile='iSGS_3364_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20230922_gpa_scRNA_hg38_test_deconx/decontX/result/iSGS_3364_1')

### Parameter setting end ###

source("scRNA_func.r")
library("SingleCellExperiment")
library("celda")
library("ggplot2")
library(patchwork)

sample_map = read_file_map("fileList1.txt")
raw_map = read_file_map("fileList2.txt")

myoptions = read_file_map("fileList3.txt", do_unlist=FALSE)
myoptions$remove_decontX = is_one(myoptions$remove_decontX)
myoptions$remove_decontX_by_contamination = as.numeric(myoptions$remove_decontX_by_contamination)

read_sce<-function(countfile){
  cat("  read", countfile, "\n")
  lst = read_scrna_data(countfile)
  counts<-lst$counts
  sce <- SingleCellExperiment(assays=list(counts=counts))
  return(sce)
}

draw_umap<-function(sce, png_file){
  df = data.frame(sce@metadata$decontX$estimates$all_cells$UMAP)
  df$contamination = unlist(sce@metadata$decontX$estimates$all_cells$contamination)

  #if we subset sce, metadata would not be changed, so we need to subset df as well
  df = df[colnames(sce),]

  g1<-ggplot(df, aes(DecontX_UMAP_1, DecontX_UMAP_2)) + geom_point(aes(color=contamination), size=1) + scale_color_gradient(low="lightgray", high="red", limits = c(0,1)) + theme_bw3() + theme(aspect.ratio = 1)
  g2<-ggplot(df, aes(x = contamination)) + geom_histogram(aes(y = after_stat(density)), color="black", fill=NA, bins = 50) + geom_density() + theme_bw3() + xlim(0,1)
  g<-g1+g2+plot_layout(ncol=2)
  png(png_file, width=2300, height=1000, res=300)
  print(g)
  dev.off()
}

sample_name=names(sample_map)[1]
for (sample_name in names(sample_map)){
  cat(sample_name, "\n")
  countfile = sample_map[sample_name]
  raw_countfile = raw_map[sample_name]

  sce <- read_sce(countfile)
  sce.raw <- read_sce(raw_countfile)

  sce <- decontX(sce, background = sce.raw)

  draw_umap(sce, paste0(sample_name, ".decontX.png"))

  meta=colData(sce)
  saveRDS(meta, paste0(sample_name, ".decontX.meta.rds"))

  if(myoptions$remove_decontX & myoptions$remove_decontX_by_contamination > 0){
    decontX_sce = sce[, sce@metadata$decontX$estimates$all_cells$contamination < myoptions$remove_decontX_by_contamination]
  
    draw_umap(decontX_sce, paste0(sample_name, ".decontX.after.png"))

    filtered = data.frame("Pre_filter_cell" = ncol(sce), "Post_filter_cell" = ncol(decontX_sce))
    write.csv(filtered, paste0(sample_name, ".decontX.filtered.csv"), row.names=FALSE)

    counts = counts(decontX_sce)
  }else{
    counts = ceiling(decontXcounts(sce))
  }
  saveRDS(counts, paste0(sample_name, ".decontX.counts.rds"))
}
