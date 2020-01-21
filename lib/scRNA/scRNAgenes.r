
library(Seurat)
library(ggplot2)

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
  g1<-DimPlot(obj, reduction = "umap", label=T, group.by="seurat_cellactivity_clusters", cols=seurat_cellactivity_colors) + NoLegend()
  l<-list(umap=g1)
  
  sample=samples[1]
  for (sample in samples){
    cells<-colnames(obj)[obj$orig.ident==sample]
    g2<-FeaturePlot(object = obj, cells = cells, features=gene) + ggtitle(sample)
    l[[sample]]=g2
  }

  sh<-floor (sqrt(length(l)))
  sw<-ceiling(length(l) / sh)
  png(filename=paste0(outFile, ".", gene, ".png"), width=sw * 3000, height=sh * 3000, res=300)
  print(CombinePlots(l, nrow=2))
  dev.off()
}
