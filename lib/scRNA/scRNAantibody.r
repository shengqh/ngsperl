library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)

obj<-finalList$obj

clusterDf<-read.csv(parFile3, stringsAsFactors = F)
clusters<-clusterDf[,cluster_name]

caCount<-table(clusters)
clusterDf$caCount<-caCount[clusters]

clusterDf<-clusterDf[order(clusterDf$caCount, clusterDf[,cluster_name], decreasing = c(1,0)),]

clusters<-factor(clusters, levels=unique(clusterDf[,cluster_name]))

obj$final_seurat_clusters<-clusters

adt<-obj[["ADT"]]
genes<-rownames(adt@data)
ncol=ceiling(sqrt(1 + length(genes)))
nrow=ceiling((1 + length(genes)) / ncol)

p1<-MyDimPlot(obj, reduction = "umap", label=T, group.by="final_seurat_clusters") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))
lst<-list("Cluster" = p1)

gene=genes[1]
for(gene in genes){
  name<-gsub("-TotalSeqC", "", gene)
  name<-gsub("TotalSeqC-", "", name)
  pgene<-FeaturePlot(object = obj, features=paste0("adt_", gene))  + NoLegend() + ggtitle(name) + theme(plot.title = element_text(hjust=0.5))
  lst[[name]] = pgene
}

png(filename=paste0(outFile, ".antibody.png"), width= ncol * 3000, height=nrow * 3000, res=300)
g<-ggarrange(plotlist=lst, ncol=ncol, nrow=nrow)
print(g)
dev.off()

png(filename=paste0(outFile, ".dot.png"), width=max(length(genes) * 100, 5000), height=2500, res=300)
p<-DotPlot(obj, assay = assay, group.by="final_seurat_clusters", features=paste0("adt_", genes), cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
  xlab("genes")
print(p)
dev.off()
