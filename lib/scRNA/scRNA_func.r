
read_cell_cluster_file<-function(fileName, sort_cluster_name="seurat_clusters"){
  result<-read.csv(fileName, stringsAsFactors = F, row.names = 1)
  
  display_sort_cluster_name = paste0("display_", sort_cluster_name)
  result[,display_sort_cluster_name] = paste0("Cluster ", result[,sort_cluster_name])
  
  cluster_names=colnames(result)[grepl("_clusters", colnames(result))]
  
  sort_clusters_num = length(unique(result[,sort_cluster_name]))
  for(cluster_name in cluster_names){
    cluster_num = length(unique(result[,cluster_name]))
    if(cluster_name == sort_cluster_name){
      next
    }

    if (cluster_num != sort_clusters_num) {
      next
    }
      
    cf<-unique(result[,c(sort_cluster_name, cluster_name)])
    if(nrow(cf) != sort_clusters_num){
      next
    }
    
    cf<-cf[order(as.numeric(cf[,sort_cluster_name]), decreasing = T),]
    cf_levels=cf[,cluster_name]
    result[,cluster_name] = factor(result[,cluster_name], levels=cf_levels)
  }
  return(result)
}

find_markers<-function(object, by_sctransform, ident.1, ident.2, min.pct = 0.5, logfc.threshold = 0.6){
  assay=ifelse(by_sctransform, "SCT", "RNA")
  markers=FindMarkers(object, assay=assay, ident.1=ident.1, ident.2=ident.2, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  markers=markers[markers$p_val_adj < 0.01,]
  return(markers)
}

