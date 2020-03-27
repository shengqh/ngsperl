
library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
genes<-unlist(strsplit(genes, ";"))
genes<-unique(genes)

obj<-finalList$obj

missgenes<-genes[!(genes %in% rownames(obj))]
write.table(data.frame("Gene"=missgenes), file="missed_genes.txt", row.names = F, col.names=F)

genes<-genes[genes %in% rownames(obj)]

clusterDf<-read.csv(parFile3, stringsAsFactors = F)
clusters<-clusterDf[,cluster_name]

caCount<-table(clusters)
clusterDf$caCount<-caCount[clusters]

clusterDf<-clusterDf[order(clusterDf$caCount, clusterDf$seurat_clusters, decreasing = c(1,0)),]

clusters<-factor(clusters, levels=unique(clusterDf[,cluster_name]))

obj$final_seurat_clusters<-clusters

samples<-unique(obj$orig.ident)

gene=genes[1]
ncol=ceiling(sqrt(1 + length(samples)))
nrow=ceiling((1 + length(samples)) / ncol)
for(gene in genes){
  p1<-DimPlot(obj, reduction = "umap", label=T, group.by="final_seurat_clusters") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))
  lst<-list("Cluster" = p1)
  sample<-samples[1]
  for (sample in samples){
    sobj<-subset(obj, cells=colnames(obj)[obj$orig.ident==sample])
    pgene<-FeaturePlot(object = obj, features=gene)  + NoLegend() + ggtitle(sample) + theme(plot.title = element_text(hjust=0.5))
    lst[[sample]]=pgene
  }

  png(filename=paste0(outFile, ".", gene, ".png"), width= ncol * 3000, height=nrow * 3000, res=300)
  g<-ggarrange(plotlist=lst, ncol=ncol, nrow=nrow)
  print(g)
  dev.off()
}

png(filename=paste0(outFile, ".dot.png"), width=max(length(genes) * 100, 5000), height=2500, res=300)
p<-DotPlot(obj, group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
  xlab("genes")
print(p)
dev.off()
