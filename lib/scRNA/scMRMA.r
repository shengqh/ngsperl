library(Seurat)
library(scMRMA)
library(plyr)

finalList<-readRDS(parFile1)
obj<-finalList$obj
options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
options<-split(options_table$V1, options_table$V2)

species=options$species
prefix<-options$prefix
db=options$db
p=as.numeric(options$p)
combined_mode=as.numeric(options$combined_mode)

res_all=get_annotation(input =  obj,species =  species,db =  db, p = p, cluster_func = scMRMA:::get_seurat_cluster, skip_first_layer = T)
res=res_all$multiR$annotationResult
res$UniformR=as.character(res_all$uniformR$annotationResult$UniformR)

cluster_names=c()

layer=colnames(res)[6]
for(layer in colnames(res)){
  curlayer = res[[layer]]
  curcelltype=plyr::mapvalues(colnames(obj), rownames(res), curlayer)

  layer_name = paste0("scMRMA_", layer)
  obj[[layer_name]] <- curcelltype

  png(paste0(prefix, layer_name, ".png"), width=4000, height=3000, res=300)
  print(DimPlot(obj,reduction = "umap",group.by = layer_name,label = TRUE,repel = TRUE)+ggplot2::ggtitle(layer_name))
  dev.off()

  layer_cluster = paste0(layer_name, "_cluster")
  obj[[layer_cluster]] <- paste0(obj$seurat_clusters, ":", curcelltype)

  png(paste0(prefix, layer_cluster, ".png"), width=4000, height=3000, res=300)
  print(DimPlot(obj,reduction = "umap",group.by = layer_cluster,label = TRUE,repel = TRUE)+ggplot2::ggtitle(layer_name))
  dev.off()
  #write.csv(layerdata$genecount, file=paste0(prefix, layer_name, ".csv"), row.names=F)

  cluster_names=c(cluster_names, layer_name, layer_cluster)
}

finalListFile=paste0(prefix, ".scMRMA.rds")
saveRDS(finalList, file=finalListFile)

clusters<-data.frame("cell" = c(1:length(obj$seurat_clusters)), "seurat_clusters"=as.numeric(as.character(obj$seurat_clusters)), stringsAsFactors = F)
rownames(clusters)<-names(obj$seurat_clusters)
for(cluster_name in cluster_names){
  clusters[[cluster_name]] = obj[[cluster_name]]
}
write.csv(clusters, file=paste0(prefix, ".cluster.csv"))
