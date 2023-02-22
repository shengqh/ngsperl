rm(list=ls()) 
outFile='doublets'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/wanjalla_lab/projects/20230220_scRNA_P8008/P8008_CW2_compressed.rds'
parFile2='/data/wanjalla_lab/projects/20230221_doublets/seurat_edgeR_betweenCluster_byCell/result/doublets.edgeR.files.csv'
parFile3=''


setwd('/data/wanjalla_lab/projects/20230221_doublets/seurat_edgeR_betweenCluster_byCell_dotplot/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(ggpubr)

obj<-read_object(parFile1)
sigGeneList<-read.csv(parFile2)
sigGeneList$designFile<-paste0(dirname(parFile2), "/", sigGeneList$designFile)
sigGeneList$sigFile<-paste0(dirname(parFile2), "/", sigGeneList$sigFile)

if(parFile3 == ""){
  cell_df<-obj@meta.data
}else{
  cell_df<-read_cell_cluster_file(parFile3)
}

params_def=read.table(parSampleFile1)
params <- setNames(as.character(params_def$V1), params_def$V2)
gene_number=as.numeric(params['gene_number'])
cluster_name=params['cluster_name']
display_cluster_name=params['display_cluster_name']

if(display_cluster_name == "seurat_clusters"){
  display_cluster_name = "display_seurat_clusters"
}
obj[["final_seurat_clusters"]]=cell_df[,display_cluster_name]

idx=1
for(idx in c(1:nrow(sigGeneList))) {
  comparison = sigGeneList$comparison[idx]

  designFile = sigGeneList$designFile[idx]
  design=read.csv(designFile)
  cells=design$Cell

  sigFile = sigGeneList$sigFile[idx]
  sig=read.csv(sigFile)
  genes=sig$X[sig$logFC > 0]
  genes=genes[1:min(gene_number, length(genes))]

  subobj=subset(obj,cells=cells)

  clusternames=unique(subobj[["final_seurat_clusters"]])

  pdf(file=paste0(comparison, ".dot.pdf"), width=max(length(genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
  p<-DotPlot(subobj, assay = "RNA", group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab("genes") + ggtitle(comparison) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Up-regulated genes") + ylab("")
  print(p)
  dev.off()
}
