
setwd('C:/projects/scratch/cqs/shengq2/alexander_gelbard_projects/20210611_NoAAC_iSGS_3364_3800_6363_scRNA/seurat_harmony_celltype_gene_localization_map/result')

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)

obj<-finalList$obj

groups_tbl<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F)
groups=split(groups_tbl$V2, groups_tbl$V1)
obj$group = unlist(groups[obj$orig.ident])

ngroup=length(unique(groups_tbl$V2))

genes_tbl<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
gene=genes_tbl$V1[1]
for (gene in genes_tbl$V1){
  png(filename=paste0(outFile, gene, ".feature_plot.png"), width= ngroup * 2000, height=2000, res=300)
  g<-FeaturePlot(obj, gene, split.by = "group", pt.size=1)
  print(g)
  dev.off()
}
