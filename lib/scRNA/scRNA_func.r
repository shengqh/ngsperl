require(data.table)

read_cell_markers_file<-function(panglao5_file, species, remove_subtype_of="", HLA_panglao5_file=""){
  #preparing cell activity database
  marker<-data.frame(fread(panglao5_file))
  if(remove_subtype_of != ""){
    remove_subtype_of<-unlist(strsplit(remove_subtype_of, ","))
    pangdb_ct <- read.table(HLA_panglao5_file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    layer="Layer2"
    removed<-c()
    for(layer in c("Layer1","Layer2", "Layer3")){
      layer_ct<-unique(pangdb_ct[,layer])
      for(rs in remove_subtype_of){
        if(rs %in% removed){
          next
        }
        
        if(rs %in% layer_ct){
          subdb<-rownames(pangdb_ct)[pangdb_ct[,layer]==rs]
          if(rs %in% subdb){ 
            subdb<-subdb[subdb != rs]
            marker<-marker[!(marker$cell.type %in% subdb),]
          }else{
            marker$cell.type[marker$cell.type %in% subdb] = rs
          }
          removed<-c(removed, rs)
        }
      }
    }
  }
  
  hsind<-regexpr(species,marker$species)
  marker_species<-marker[hsind>0 & marker$ubiquitousness.index<0.05,]
  if (species=="Mm") {
    ##change the gene symbol only keep the first letter capitalize
    marker_species$official.gene.symbol<-toMouseGeneSymbol(marker_species$official.gene.symbol)
  }
  if (species=="Hs") {
    marker_species$official.gene.symbol<-toupper(marker_species$official.gene.symbol)
  }
  cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  return(list(cellType=cellType, weight=weight))
}

toMouseGeneSymbol<-function(x){
  result=paste0(toupper(substr(x,1,1)),tolower(substr(x,2,nchar(x))))
  return(result)
}

read_cell_markers_file<-function(panglao5_file, species, remove_subtype_of="", HLA_panglao5_file=""){
  #preparing cell activity database
  marker<-data.frame(fread(panglao5_file))
  if(remove_subtype_of != ""){
    pangdb_ct <- read.table(HLA_panglao5_file,header = F,row.names = 1,sep = "\t",stringsAsFactors = F)

  }

  hsind<-regexpr(species,marker[,1])
  marker_species<-marker[hsind>0 & marker$ubiquitousness.index<0.05,]
  if (species=="Mm") {
    ##change the gene symbol only keep the first letter capitalize
    marker_species$official.gene.symbol<-toMouseGeneSymbol(marker_species$official.gene.symbol)
  }
  if (species=="Hs") {
    marker_species$official.gene.symbol<-toupper(marker_species$official.gene.symbol)
  }
  cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  finalList$cell_activity_database<-list(cellType=cellType, weight=weight)
}

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

get_cluster_count<-function(counts, clusters){
  if(is.null(levels(clusters))){
    allClusters=unique(clusters)
  }else{
    allClusters=levels(clusters)
  }
  cluster<-allClusters[1]
  csums=lapply(allClusters, function(cluster){
    #cat(cluster, "\n")
    cells=names(clusters)[clusters==cluster]
    subcounts=counts[,cells]
    #cat(ncol(subcounts), "\n")
    Matrix::rowSums(subcounts)
  })
  
  result=do.call(cbind, csums)
  colnames(result)=allClusters
  
  gcount=Matrix::rowSums(result)
  result=result[gcount > 0,]
  return(result)
}

get_group_count=function(curobj, groupName="active.ident") {
  counts=GetAssayData(curobj,assay="RNA",slot="counts")
  curgroups=curobj[[groupName]]
  clusters=curgroups[,1]
  names(clusters)=rownames(curgroups)
  result=get_cluster_count(counts, clusters)
  return(result)
}
