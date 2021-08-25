library(Seurat)
library(ggplot2)
library(ggpubr)

cell_df<-read_cell_cluster_file(parFile3)

finalList<-readRDS(parFile1)
sigGeneList<-read.csv(parFile2)
sigGeneList$designFile<-paste0(dirname(parFile2), "/", sigGeneList$designFile)
sigGeneList$sigFile<-paste0(dirname(parFile2), "/", sigGeneList$sigFile)

obj<-finalList$obj
assay=ifelse("SCT" %in% names(obj@assays), "SCT", "RNA")

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
  p<-DotPlot(subobj, assay = assay, group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab("genes") + ggtitle(comparison) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Up-regulated genes") + ylab("")
  print(p)
  dev.off()
}
