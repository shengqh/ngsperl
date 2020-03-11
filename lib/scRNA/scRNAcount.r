
library(Seurat)

finalList<-readRDS(parFile1)

output_cluster<-function(obj, outFile, cluster_name){
  clusterCells<-obj[[cluster_name]]
  cluster<-data.frame(Cluster=unique(clusterCells[,1]))
  files<-gsub(" ", "_", cluster$Cluster)
  files<-gsub(":", "_", files)
  files<-gsub("_+", "_", files)
  cluster$DE<-files

  write.csv(cluster, file=paste0(outFile, ".count.files.csv"), row.names=F, quote=F)
  
  des_unique<-unique(cluster$DE)
  de<-des_unique[1]
  for(de in des_unique){
    de_ids<-cluster$Cluster[cluster$DE == de]
    cells<-rownames(clusterCells)[clusterCells[,1] %in% de_ids]
    de_obj<-obj[, cells]
    de_count<-as.matrix(de_obj[["RNA"]]@counts)
    de_file_name<-paste0(outFile, ".", de, ".count")
    saveRDS(de_count, paste0(de_file_name,".rds"))
    
    orig.ident<-de_obj[["orig.ident"]]
    colnames(orig.ident)<-c("Sample")
    write.csv(orig.ident, paste0(de_file_name,".sample.csv"), quote=F)
  }
}

obj<-finalList$obj

clusterDf<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
obj[[cluster_name]]<-clusterDf[names(obj$orig.ident), cluster_name]

output_cluster(obj, outFile, cluster_name)
