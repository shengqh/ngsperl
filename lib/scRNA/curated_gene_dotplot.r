
library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
obj<-finalList$obj

assay=ifelse("SCT" %in% names(obj@assays), "SCT", "RNA")

genes_df=read.table(parSampleFile1, sep="\t", stringsAsFactors=F)

cluster_df=read.table(parSampleFile2, sep="\t", stringsAsFactors=F)
cluster_request <- tapply(cluster_df$V1,cluster_df$V2,list)

params_def=read.table(parSampleFile3, sep="\t", stringsAsFactors=F)
params=tapply(params_def$V1,params_def$V2,list)

cell_df<-read_cell_cluster_file(parFile2)

obj[["final_seurat_clusters"]]=cell_df[,params$display_cluster_name]

#assay=ifelse(params$by_sctransform=="1", "SCT", "RNA")
assaydata=GetAssayData(obj, assay=assay)
allgenes=rownames(assaydata)
rm(assaydata)
genes_df=subset(genes_df, genes_df$V1 %in% allgenes)

is_one_cluster = length(cluster_request) == 1
idx=1
for(idx in c(1:length(cluster_request))) {
  geneset_name=names(cluster_request)[idx]
  clusternames = as.character(cluster_request[[idx]])

  if(clusternames=='all'){
    subobj=obj
  }else{
    cells=rownames(cell_df)[cell_df[,cluster_name] %in% clusternames]
    subobj=subset(obj, cells=cells)
  }

  cname_df=data.frame("cluster_name"=subobj[[params$cluster_name]], "display_cluster_name"=subobj[["final_seurat_clusters"]])

  cur_genes_df=subset(genes_df, genes_df$V2 == geneset_name)
  genes=unique(cur_genes_df$V1)

  pdf(file=paste0(geneset_name, ".dot.pdf"), width=max(length(genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
  p<-DotPlot(subobj, assay = assay, group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    theme(plot.title = element_text(hjust = 0.5)) + xlab(gsub("_", " ", geneset_name)) + ylab("")
  
  if(!is_one_cluster){
    p<-p+ggtitle(geneset_name)
  }
  
  print(p)
  dev.off()
}
