library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
obj<-finalList$obj

genes_df=read.table(parSampleFile1, sep="\t")

cluster_df=read.table(parSampleFile2, sep="\t")
cluster_request <- tapply(cluster_df$V1,cluster_df$V2,list)

params_def=read.table(parSampleFile3, sep="\t")
params <- setNames(as.character(params_def$V1), params_def$V2)
gene_number=as.numeric(params['gene_number'])
cluster_name=params['cluster_name']
display_cluster_name=params['display_cluster_name']

cell_df<-read_cell_cluster_file(parFile2)

obj[["final_seurat_clusters"]]=cell_df[,display_cluster_name]

is_one_cluster = length(cluster_request) == 1
idx=1
for(idx in c(1:length(cluster_request))) {
  curname=names(cluster_request)[idx]
  clusternames = as.character(cluster_request[[idx]])

  cells=rownames(cell_df)[cell_df[,cluster_name] %in% clusternames]
  subobj=subset(obj, cells=cells)

  cname_df=data.frame("cluster_name"=subobj[[cluster_name]], "display_cluster_name"=subobj[["final_seurat_clusters"]])

  cur_genes_df=subset(genes_df, genes_df$V3 == curname)
  dms = tapply(cur_genes_df$V1, cur_genes_df$V2, list)

  idy=1
  for(idy in c(1:length(dms))) {
    geneset_name = names(dms)[idy]
    genes=unique(unname(unlist(dms[idy])))
    genes=genes[order(genes)]
    
    fgenes<-FetchData(subobj, genes)
    valid_genes<-genes[genes %in% colnames(fgenes)]
    
    pdf(file=paste0(curname, ".", gsub(" ", "_", geneset_name), ".dot.pdf"), width=max(length(valid_genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
    p<-DotPlot(subobj, group.by="final_seurat_clusters", features=valid_genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
      theme(plot.title = element_text(hjust = 0.5)) + xlab(gsub("_", " ", geneset_name)) + ylab("")
    
    if(!is_one_cluster){
      p<-p+ggtitle(curname)
    }
    
    print(p)
    dev.off()
  }
}
