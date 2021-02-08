library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
sigGeneList<-read.csv(parFile2)
sigGeneList$designFile<-paste0(dirname(parFile2), "/", sigGeneList$designFile)
sigGeneList$sigFile<-paste0(dirname(parFile2), "/", sigGeneList$sigFile)

obj<-finalList$obj

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

  subobj[["final_seurat_clusters"]]=subobj[[cluster_name]]

  clusternames=unique(subobj[["final_seurat_clusters"]])

  pdf(file=paste0(comparison, ".dot.pdf"), width=max(length(genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
  p<-DotPlot(subobj, assay = "RNA", group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab("genes") + ggtitle(comparison) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Up-regulated genes") + ylab(cluster_name)
  print(p)
  dev.off()
}
