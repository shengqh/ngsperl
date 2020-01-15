library(Seurat)

finalList<-readRDS(parFile1)

output_cluster<-function(obj, cts, cts_name, DE_by_celltype=TRUE){
  cts_cluster<-data.frame(Cluster=levels(obj), CellType=names(cts), Score=cts )
  cts_files<-gsub(" ", "_", cts_cluster$CellType)
  if(DE_by_celltype){
    cts_cluster$DE<-cts_files
  }else{
    cts_cluster$DE<-paste0(cts_cluster$Cluster, "_", cts_files)
  }  
  
  write.csv(cts_cluster, file=paste0(cts_name, ".cluster.csv"), row.names=F, quote=F)
  
  des_unique<-unique(cts_cluster$DE)
  for(de in des_unique){
    de_ids<-cts_cluster$Cluster[cts_cluster$DE == de]
    de_obj<-subset(obj, ident=de_ids)
    de_count<-as.matrix(de_obj[["RNA"]]@counts)
    de_file_name<-paste0(cts_name, ".", de, ".count")
    saveRDS(de_count, paste0(de_file_name,".rds"))
    
    orig.ident<-de_obj[["orig.ident"]]
    colnames(orig.ident)<-c("Sample")
    write.csv(orig.ident, paste0(de_file_name,".sample.csv"), quote=F)
  }
}

cts<-finalList$cell_activity_database$predicted$max_cta
obj<-finalList$obj
output_cluster(obj, cts, outFile, DE_by_celltype)
