
library(Seurat)
library(ggplot2)

finalList<-readRDS(parFile1)
geneFile<-parFile2

obj<-finalList$obj
seurat_colors<-finalList$seurat_colors
seurat_cellactivity_colors<-finalList$seurat_cellactivity_colors

celltypes<-read.table(geneFile, sep="\t", header=T, stringsAsFactors = F)
celltypes$Gene<-gsub("\\s.+", "", celltypes$Gene)
celltypes$Gene<-toupper(celltypes$Gene)

missGenes<-celltypes[!(celltypes$Gene %in% rownames(obj)),]
write.csv(missGenes, "miss_gene.csv")
celltypes<-celltypes[celltypes$Gene %in% rownames(obj),]

#allgenes<-data.frame(Gene=rownames(obj))
#write.csv(allgenes, "all_gene.csv")
#celltypes$Gene<-paste0(substr(celltypes$Gene,1,1),substr(tolower(celltypes$Gene),2,nchar(celltypes$Gene)))

max_cta<-finalList$cell_activity_database$predicted$max_cta
cta_df<-data.frame(Cluster=c(1:length(max_cta))-1, Celltype=names(max_cta), stringsAsFactors = F)
ct_unique<-unique(cta_df$Celltype)

ct<-unique(celltypes$Subtype)[1]
for(ct in unique(celltypes$Subtype)){
  subcelltypes<-celltypes[celltypes$Subtype==ct,]
  expressedGenes<-subcelltypes$Gene[subcelltypes$Status=="expressed"]
  absentGenes<-toupper(subcelltypes$Gene[subcelltypes$Status=="absent"])
  
  png(filename=paste0(outFile, ".", ct, ".png"), width=5000, height=2500, res=300)
  p<-DotPlot(obj, group.by="seurat_cellactivity_clusters", features=expressedGenes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab(paste0(ct, " expressed genes"))
  if(length(absentGenes) > 0){
    e<-DotPlot(obj, group.by="seurat_cellactivity_clusters", features=absentGenes, cols = c("lightgrey", "blue"), dot.scale = 8) + RotatedAxis()+
      xlab(paste0(ct, " absent genes"))
    p2<-CombinePlots(plots = list(p, e))
    print(p2)
  }else{
    print(p)
  }
  dev.off()
}
