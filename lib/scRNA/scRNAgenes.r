
library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
genes<-unlist(strsplit(genes, ";"))

obj<-finalList$obj
seurat_colors<-finalList$seurat_colors
seurat_cellactivity_colors<-finalList$seurat_cellactivity_colors

samples<-unique(obj$orig.ident)

missgenes<-genes[!(genes %in% rownames(obj))]
write.table(data.frame("Gene"=missgenes), file="missed_genes.txt", row.names = F, col.names=F)

genes<-genes[genes %in% rownames(obj)]

gene=genes[1]
for(gene in genes){
  p1<-DimPlot(obj, reduction = "umap", label=T, group.by="seurat_cellactivity_clusters", cols=seurat_cellactivity_colors) + NoLegend()
  p2<-FeaturePlot(object = obj, features=gene, split.by = "orig.ident")

  png(filename=paste0(outFile, ".", gene, ".png"), width= (length(samples) + 1) * 3000, height=3000, res=300)
  g<-ggarrange(p1, p2, widths = c(1,length(samples)))
  print(g)
  dev.off()
}
