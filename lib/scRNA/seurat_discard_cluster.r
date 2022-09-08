library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)

options(future.globals.maxSize= 10779361280)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
reduction<-myoptions$reduction
assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")

dynamic_discard_cluster=unlist(strsplit(myoptions$dynamic_discard_cluster, ','))

cluster_layer=myoptions$cluster_layer
celltype_layer=myoptions$celltype_layer
cluster_celltype_layer=myoptions$cluster_celltype_layer

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

prefix<-outFile

if(!exists("obj")){
  obj=read_object(parFile1, parFile2)
  Idents(obj)<-cluster_celltype_layer
}

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  bubble_genes<-unique(genes_df$gene)
}

output_figure<-function(obj, cluster_layer, celltype_layer, cluster_celltype_layer, has_bubblemap, bubblemap_file, file_name){
  meta.data = obj@meta.data
  ct<-meta.data[!duplicated(meta.data[,cluster_layer]),]
  ct<-ct[order(ct[,cluster_layer]),]
  
  g<-DimPlot(obj, group.by = cluster_layer, label=T) + ggtitle(celltype_layer)+
    scale_color_discrete(labels = ct[,cluster_celltype_layer])
  if(has_bubblemap){
    g<-g+get_bubble_plot(obj, cluster_layer, celltype_layer, bubblemap_file, assay="RNA")
    g<-g+plot_layout(ncol = 2, widths = c(4, 6))
    width=11000
  }else{
    width=4300
  }
  
  png(file_name, width=width, height=4000, res=300)
  print(g)
  dev.off()
}

output_figure(obj, cluster_layer, celltype_layer, cluster_celltype_layer, has_bubblemap, bubblemap_file, paste0(prefix, ".umap.old.png") )

meta.data=obj@meta.data
meta.data=meta.data[!(meta.data[,cluster_layer] %in% dynamic_discard_cluster),]
obj=subset(obj, cells=rownames(meta.data))
obj@meta.data = meta.data

output_figure(obj, cluster_layer, celltype_layer, cluster_celltype_layer, has_bubblemap, bubblemap_file, paste0(prefix, ".umap.png") )

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))
