rm(list=ls())
sample_name='S01_ClassPTC_BRAF'
outFile='S01_ClassPTC_BRAF'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/12904_RB_VisiumHD/20251014_12904_VisiumHD_cellsegment/SingleR/result/S01_ClassPTC_BRAF')

### Parameter setting end ###

library(sf)
library(SingleR)
library(Seurat)
library(ggplot2)
library(patchwork)
library(celldex)
library(logger)

source("scRNA_func.r")

log_appender(appender_tee(paste0(sample_name, ".log")))

cp = as.character(packageVersion("celldex"))

major_minor <- function(x, split = "."){
  s <- strsplit(x, split = split, fixed = TRUE)
  y <- do.call(rbind.data.frame, s)
  setNames(y, c("major", "minor"))
}
cpv=major_minor(cp)
if(cpv[1,1] == "1" && as.numeric(cpv[1,2]) < 13){
  stop("Please install celldex version 1.13.0 or later")
}

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

bubblemap_width=to_numeric(myoptions$bubblemap_width, 4000)
bubblemap_height=to_numeric(myoptions$bubblemap_height, 2000)
bubblemap_unit=ifelse(bubblemap_width > 50, "px", "in")

assay=myoptions$assay
is_polygons=assay == "Spatial.Polygons"

if(myoptions$species == "Mm"){
  ct_ref = MouseRNAseqData()
}else if (myoptions$species == "Hs") {
  ct_ref <- HumanPrimaryCellAtlasData()
}else{
  log_warn(paste0("Cannot find singleR ref db for species ", myoptions$species, ", use HumanPrimaryCellAtlasData"))
  ct_ref <- HumanPrimaryCellAtlasData()
}

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

if(is_polygons){
  min_umi=as.numeric(myoptions$nCount_cutoff)
  log_info(paste0("Subsetting polygons with min UMIs ", min_umi, " ..."))
  DefaultAssay(obj) <- "Spatial.Polygons"
  obj <- subset(obj, subset = nCount_Spatial.Polygons >= min_umi)
}

force=TRUE
rds_file=paste0(outFile, ".SingleR.rds")
if(file.exists(rds_file) & !force){
  labels<-readRDS(rds_file)
}else{
  log_info("Converting to SingleCellExperiment object...")
  sce=as.SingleCellExperiment(DietSeurat(obj))
  log_info("Running SingleR...")
  labels<-SingleR(sce, ref=ct_ref, assay.type.test=1, labels=ct_ref$label.main)
  rm(sce)

  labels$pruned.labels[is.na(labels$pruned.labels)] <- "unclassified"
  saveRDS(labels, file=rds_file)
}

stopifnot(all(colnames(obj) == rownames(labels)))

ct_tbl<-table(labels$pruned.labels)
ct_tbl<-ct_tbl[names(ct_tbl) != "unclassified"]
ct_tbl<-ct_tbl[order(ct_tbl, decreasing = T)]
ct_tbl<-ct_tbl[ct_tbl >= sum(ct_tbl) * 0.01]

major_cts=names(ct_tbl)

labels$major_labels=labels$pruned.labels
labels$major_labels[!(labels$major_labels %in% c(major_cts, "unclassified"))]="other"

obj <- AddMetaData(obj, metadata = labels$major_labels, col.name = "SingleR_major_labels")

ct_name="SingleR_labels"
obj <- AddMetaData(obj, metadata = labels$pruned.labels, col.name = ct_name)

log_info("Saving meta data...")
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

df<-data.frame("SingleR"=obj$SingleR_labels, "Sample"=obj$orig.ident)
df_tbl<-table(df$SingleR,df$Sample)
write.csv(df_tbl, paste0(outFile, ".SingleR_Sample.csv"))

major_obj<-subset(obj, cells=colnames(obj)[obj$SingleR_major_labels %in% major_cts])
rm(obj)

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

major_obj=get_category_with_min_percentage(major_obj, ct_name, 0.01)
ct_name_count = paste0(ct_name, "_count")
major_obj@meta.data = add_column_count(major_obj@meta.data, ct_name, ct_name_count)

log_info("Visualizing...")
if("umap" %in% names(major_obj@reductions)){
  g=get_dim_plot_labelby(major_obj, label.by = ct_name, reduction="umap", pt.size=0.1) + theme(plot.title=element_blank())
  ggsave(paste0(outFile, ".SingleR.png"), g, width=6, height=4, units="in", dpi=300, bg="white")
}

if(has_bubblemap){
  if(is_polygons){
    assay="Spatial.Polygons"
  }else{
    assay="RNA"
  }
  g<-get_bubble_plot(
    obj=major_obj, 
    cur_res=NA, 
    cur_celltype=ct_name_count, 
    bubblemap_file, 
    assay=assay, 
    species=myoptions$species,
    dot.scale=4)
  ggsave(paste0(outFile, ".SingleR.dot.png"), g, width=bubblemap_width, height=bubblemap_height, units=bubblemap_unit, dpi=300, bg="white")
}
rm(major_obj)

slim_labels=subset(labels, labels %in% major_cts)
nct=length(unique(slim_labels$pruned.labels))
ncol=ceiling(sqrt(nct))
nrow=ceiling(nct / ncol)

png(paste0(outFile, ".delta.png"), width=max(1500, ncol * 500), height=max(1000, nrow * 500), res=300)
plotDeltaDistribution(slim_labels, ncol = ncol)
dev.off()

png(paste0(outFile, ".score.png"), width=4000, height=3000, res=300)
plotScoreHeatmap(slim_labels)
dev.off()

log_info("Done.\n")
