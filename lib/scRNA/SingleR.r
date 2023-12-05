rm(list=ls()) 
sample_name='DMSO_F_1513'
outFile='DMSO_F_1513'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20231130_scRNA_8870_mouse_redo_cellbender/raw_qc_sct2_SingleR/result/DMSO_F_1513')

### Parameter setting end ###

source("scRNA_func.r")
library(SingleR)
library(Seurat)
library(ggplot2)
library(patchwork)
library(celldex)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

if(myoptions$species == "Mm"){
  ct_ref = MouseRNAseqData()
}else if (myoptions$species == "Hs") {
  ct_ref <- HumanPrimaryCellAtlasData()
}else{
  warning(paste0("Cannot find singleR ref db for species ", myoptions$species, ", use HumanPrimaryCellAtlasData"))
  ct_ref <- HumanPrimaryCellAtlasData()
}

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

force=TRUE
rds_file=paste0(outFile, ".SingleR.rds")
if(file.exists(rds_file) & !force){
  labels<-readRDS(rds_file)
}else{
  sce=as.SingleCellExperiment(DietSeurat(obj))
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

ct_name="SingleR_labels"
obj <- AddMetaData(obj, metadata = labels$pruned.labels, col.name = ct_name)
ct_name="SingleR_major_labels"
obj <- AddMetaData(obj, metadata = labels$major_labels, col.name = ct_name)

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

df<-data.frame("SingleR"=obj$SingleR_labels, "Sample"=obj$orig.ident)
df_tbl<-table(df$SingleR,df$Sample)
write.csv(df_tbl, paste0(outFile, ".SingleR_Sample.csv"))

major_obj<-subset(obj, cells=colnames(obj)[obj$SingleR_major_labels %in% major_cts])
rm(obj)

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

g1=MyDimPlot(major_obj, group.by = ct_name, reduction="umap", label=T)

if(has_bubblemap){
  g2<-get_bubble_plot(major_obj, NA, ct_name, bubblemap_file, assay="RNA", 
    species=myoptions$species)
  layout <- "ABB"
  g<-g1+g2+plot_layout(design=layout)
  width=6300
}else{
  g<-g1
  width=2300
}
height=2000

ggsave(paste0(outFile, ".SingleR.png"), g, width=width, height=height, units="px", dpi=300, bg="white")
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

