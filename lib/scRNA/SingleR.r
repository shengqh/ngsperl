rm(list=ls()) 
outFile='PH_combine'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/scratch/cqs/shengq2/paula_hurley_projects/20220824_scRNA_7467_benign_hg38/seurat_sct_harmony/result/PH_combine.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/scratch/cqs/shengq2/paula_hurley_projects/20220824_scRNA_7467_benign_hg38/seurat_sct_harmony_SingleR/result')

### Parameter setting end ###

source("scRNA_func.r")
library(SingleR)
library(Seurat)
library(ggplot2)
library(patchwork)
library(celldex)

hpca.se <- HumanPrimaryCellAtlasData()

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

if(!exists("obj")){
  obj=read_object(parFile1)
}

sce=as.SingleCellExperiment(DietSeurat(obj))
labels<-SingleR(sce, ref=hpca.se, assay.type.test=1, labels=hpca.se$label.main)

saveRDS(labels, file=paste0(outFile, ".SingleR.rds"))

stopifnot(all(colnames(obj) == rownames(labels)))

ct_tbl<-table(labels$labels)
ct_tbl<-ct_tbl[order(ct_tbl, decreasing = T)]
ct_tbl<-ct_tbl[ct_tbl >= sum(ct_tbl) * 0.01]

labels$major_labels=labels$labels
labels$major_labels[!(labels$major_labels %in% names(ct_tbl))]="other"

ct_name="SingleR_labels"
obj <- AddMetaData(obj, metadata = labels$labels, col.name = ct_name)
ct_name="SingleR_major_labels"
obj <- AddMetaData(obj, metadata = labels$major_labels, col.name = ct_name)

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

ct<-data.frame("SingleR"=obj$SingleR_labels, "Sample"=obj$orig.ident)
ct_tbl<-table(ct$SingleR,ct$Sample)
write.csv(ct_tbl, paste0(outFile, ".SingleR_Sample.csv"))

major_obj<-subset(obj, cells=colnames(obj)[obj$SingleR_major_labels != "other"])

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

g1=DimPlot(major_obj, group.by = ct_name, reduction="umap", label=T)

if(has_bubblemap){
  g2<-get_bubble_plot(major_obj, NA, ct_name, bubblemap_file, assay="RNA")
  layout <- "ABB"
  g<-g1+g2+plot_layout(design=layout)
  width=6300
}else{
  g<-g1
  width=2300
}
height=2000

png(paste0(outFile, ".SingleR.png"), width=width, height=height, res=300)
print(g)
dev.off()
