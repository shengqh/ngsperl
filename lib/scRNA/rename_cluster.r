
source("scRNA_func.r")

library(Seurat)
library(ggplot2)

myoptions<-read.table(parSampleFile1, stringsAsFactors = F, sep="\t", header=F)
myoptions<-split(myoptions$V1, myoptions$V2)

newnames<-read.table(parSampleFile2, stringsAsFactors = F, sep="\t", header=F)
newnames<-split(newnames$V1, newnames$V2)

clusters<-read.csv(parFile3, stringsAsFactors = F, header=T)
rownames(clusters)=clusters$seurat_clusters

renamed_column = paste0("renamed_", myoptions$celltype_name)
clusters[,renamed_column] = clusters[,myoptions$celltype_name]
clusters[names(newnames),renamed_column] = unlist(newnames)

seurat_renamed_column=paste0("seurat_", renamed_column)
clusters[,seurat_renamed_column] = paste0(clusters$seurat_clusters, " : ", clusters[,renamed_column])
clusters[,seurat_renamed_column] = factor(clusters[,seurat_renamed_column], levels=clusters[,seurat_renamed_column])

write.csv(clusters, file=paste0(outFile, ".rename_celltype.csv"), row.names=F)

cells<-read.csv(parFile2, stringsAsFactors = F, row.names=1, header=T)

renames<-split(clusters[,renamed_column], clusters$seurat_clusters)
cells[,renamed_column]=unlist(renames[as.character(cells$seurat_clusters)])

renames<-split(clusters[,seurat_renamed_column], clusters$seurat_clusters)
cells[,seurat_renamed_column]=unlist(renames[as.character(cells$seurat_clusters)])

write.csv(cells, file=paste0(outFile, ".rename_cluster.csv"))

finalList<-readRDS(parFile1)
obj<-finalList$obj

#make sure with same cell order
cells<-cells[colnames(obj),]

obj[[seurat_renamed_column]]<-cells[,seurat_renamed_column]

png(file=paste0(outFile, ".rename_cluster.png"), width=4000, height=3000, res=300)
g<-MyDimPlot(object = obj, reduction = 'umap', label=TRUE, group.by=seurat_renamed_column) + guides(colour = guide_legend(override.aes = list(size = 3), ncol=1))
print(g)
dev.off()

cellcounts<-data.frame(Sample=obj$orig.ident, Cluster=obj[[seurat_renamed_column]])
ctable<-t(table(cellcounts))

ctable_perThousand <- round(ctable / colSums(ctable) * 1000)
colnames(ctable_perThousand)<-paste0(colnames(ctable), "_perThousand")

final<-cbind(ctable, ctable_perThousand)
write.csv(final, file=paste0(outFile, ".rename_cluster.summery.csv"))
