
library(Seurat)
library(ggplot2)

finalList<-readRDS(parFile1)
geneFile<-parFile2

obj<-finalList$obj

celltypes<-read.table(geneFile, sep="\t", header=T, stringsAsFactors = F)
celltypes$Gene<-gsub("\\s.+", "", celltypes$Gene)
celltypes$Gene<-toupper(celltypes$Gene)

missGenes<-celltypes[!(celltypes$Gene %in% rownames(obj)),]
write.csv(missGenes, "miss_gene.csv")
celltypes<-celltypes[celltypes$Gene %in% rownames(obj),]

#allgenes<-data.frame(Gene=rownames(obj))
#write.csv(allgenes, "all_gene.csv")
#celltypes$Gene<-paste0(substr(celltypes$Gene,1,1),substr(tolower(celltypes$Gene),2,nchar(celltypes$Gene)))

clusterDf<-read.csv(parFile3, stringsAsFactors = F, row.names=1)
clusterDf<-clusterDf[colnames(obj),]
clusterDf$sample<-obj$orig.ident

clusters<-clusterDf[,cluster_name]

caCount<-table(clusters)
clusterDf$caCount<-caCount[clusters]

clusterDf<-clusterDf[order(-clusterDf$caCount, clusterDf[,cluster_name]),]

clusters<-factor(clusters, levels=unique(clusterDf[,cluster_name]))

obj$final_seurat_clusters<-clusters

if (exists("samples")){
  sample_idents = unlist(strsplit(samples, ";"))
  sample_cells<-rownames(clusterDf)[clusterDf$sample %in% sample_idents]
  obj=subset(obj, cells=sample_cells)
}

ct<-unique(celltypes$Celltype)[1]
for(ct in unique(celltypes$Celltype)){
  subcelltypes<-celltypes[celltypes$Celltype==ct,]
  expressedGenes<-toupper(unique(subcelltypes$Gene[subcelltypes$Status=="expressed"]))
  absentGenes<-toupper(unique(subcelltypes$Gene[subcelltypes$Status=="absent"]))
  
  width=max(5000, nrow(subcelltypes) * 70)
  png(filename=paste0(outFile, ".", ct, ".png"), width=width, height=2500, res=300)
  p<-DotPlot(obj, group.by="final_seurat_clusters", features=expressedGenes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab(paste0(ct, " expressed genes"))
  if(length(absentGenes) > 0){
    e<-DotPlot(obj, group.by="final_seurat_clusters", features=absentGenes, cols = c("lightgrey", "blue"), dot.scale = 8) + RotatedAxis()+
      xlab(paste0(ct, " absent genes"))
    p2<-CombinePlots(plots = list(p, e))
    print(p2)
  }else{
    print(p)
  }
  dev.off()
}
