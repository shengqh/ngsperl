
library(Seurat)
library(ggplot2)

finalList<-readRDS(parFile1)
obj<-finalList$obj

newnames<-read.table(parSampleFile2, stringsAsFactors = F, sep="\t", header=F)
clusters<-data.frame("cell" = c(1:length(obj$seurat_clusters)), "seurat_clusters"=as.numeric(as.character(obj$seurat_clusters)), "cellactivity_clusters"=obj$cellactivity_clusters, "seurat_cellactivity_clusters"=obj$seurat_cellactivity_clusters, stringsAsFactors = F)
clusters<-clusters[order(clusters$seurat_clusters),]
clusters$renamed_cellactivity_clusters<-clusters$cellactivity_clusters

idx<-1
for(idx in c(1:nrow(newnames))){
  oldcluster<-newnames$V2[idx]
  newcluster<-newnames$V1[idx]
  clusters$renamed_cellactivity_clusters[clusters$seurat_clusters==oldcluster]<-newcluster
}

clusters$renamed_seurat_cellactivity_clusters<-paste0(clusters$seurat_clusters, " : ", clusters$renamed_cellactivity_clusters)
clusters$renamed_seurat_cellactivity_clusters<-factor(clusters$renamed_seurat_cellactivity_clusters, levels=unique(clusters$renamed_seurat_cellactivity_clusters))
clusters<-clusters[order(clusters$cell),]
clusters$final_clusters<-clusters$renamed_cellactivity_clusters
clusters$final_seurat_clusters<-clusters$renamed_seurat_cellactivity_clusters

write.csv(clusters, file=paste0(outFile, ".rename_cluster.csv"))

obj$renamed_cellactivity_clusters<-clusters$renamed_cellactivity_clusters
obj$renamed_seurat_cellactivity_clusters<-clusters$renamed_seurat_cellactivity_clusters

png(file=paste0(outFile, ".rename_cluster.png"), width=3200, height=3000, res=300)
g<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="renamed_seurat_cellactivity_clusters")
print(g)
dev.off()

cellcounts<-data.frame(Sample=obj$orig.ident, Cluster=obj$renamed_seurat_cellactivity_clusters)
ctable<-t(table(cellcounts))

ctable_perThousand <- round(ctable / colSums(ctable) * 1000)
colnames(ctable_perThousand)<-paste0(colnames(ctable), "_perThousand")

final<-cbind(ctable, ctable_perThousand)
write.csv(final, file=paste0(outFile, ".rename_clusters.cells.csv"))
