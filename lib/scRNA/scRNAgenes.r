
library(Seurat)
library(ggplot2)

finalList<-readRDS(parFile1)
genes<-unlist(strsplit(genes, ";"))

obj<-finalList$obj
seurat_colors<-finalList$seurat_colors
seurat_cellactivity_colors<-finalList$seurat_cellactivity_colors

for(gene in genes){
  png(filename=paste0(outFile, ".", gene, ".png"), width=5000, height=2500, res=300)
  g1<-DimPlot(obj, reduction = "umap", label=T, group.by="seurat_cellactivity_clusters", cols=seurat_cellactivity_colors) + NoLegend()
  g2<-FeaturePlot(object = obj, features=genes) 
  print(CombinePlots(list(g1, g2)))
  dev.off()
}
