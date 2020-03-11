
library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
genes<-unlist(strsplit(genes, ";"))

obj<-finalList$obj

missgenes<-genes[!(genes %in% rownames(obj))]
write.table(data.frame("Gene"=missgenes), file="missed_genes.txt", row.names = F, col.names=F)

genes<-genes[genes %in% rownames(obj)]

clusterDf<-read.csv(parFile3, stringsAsFactors = F)
clusters<-clusterDf[,cluster_name]

caCount<-table(clusters)
clusterDf$caCount<-caCount[clusters]

clusterDf<-clusterDf[order(-clusterDf$caCount, clusterDf$seurat_cluters),]

clusters<-factor(clusters, levels=unique(clusterDf[,cluster_name]))

obj$final_seurat_clusters<-clusters

samples<-unique(obj$orig.ident)

gene=genes[1]
#ncol=ceiling(sqrt(1 + length(samples)))
#nrow=ceiling((1 + length(samples)) / ncol)
for(gene in genes){
  p1<-DimPlot(obj, reduction = "umap", label=T, group.by="final_seurat_clusters") + NoLegend()
  p2<-FeaturePlot(object = obj, features=gene, split.by = "orig.ident")

  png(filename=paste0(outFile, ".", gene, ".png"), width= (length(samples) + 1) * 3000, height=3000, res=300)
  g<-ggarrange(p1, p2, widths = c(1,length(samples)))
  print(g)
  dev.off()
}

png(filename=paste0(outFile, ".dot.png"), width=5000, height=2500, res=300)
p<-DotPlot(obj, group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
  xlab("genes")
print(p)
dev.off()
