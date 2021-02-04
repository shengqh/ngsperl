library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
obj<-finalList$obj

cell_df<-read.csv(parFile2, stringsAsFactors = F, row.names = 1)
cell_df<-cell_df[colnames(obj),]

cluster_df=read.table(parSampleFile1)
cluster_request <- tapply(cluster_df$V1,cluster_df$V2,list)

params_def=read.table(parSampleFile2)
params <- setNames(as.character(params_def$V1), params_def$V2)
gene_number=as.numeric(params['gene_number'])
cluster_name=params['cluster_name']
display_cluster_name=params['display_cluster_name']

cell_df[,cluster_name]=as.character(cell_df[,cluster_name])
obj[["final_seurat_clusters"]]=cell_df[,display_cluster_name]

cnames=unique(cell_df[,c(cluster_name, display_cluster_name)])
rownames(cnames)=cnames[,cluster_name]
cnames[,display_cluster_name] = gsub('[ :]+', '_', cnames[,display_cluster_name])

idx=1
for(idx in c(1:length(cluster_request))) {
  curname=names(cluster_request)[idx]
  clusternames = as.character(cluster_request[[idx]])

  cells=rownames(cell_df)[cell_df[,cluster_name] %in% clusternames]
  subobj=subset(obj, cells=cells)

  cname_df=data.frame("cluster_name"=subobj[[cluster_name]], "display_cluster_name"=subobj[["final_seurat_clusters"]])

  markers=FindAllMarkers(subobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  write.csv(markers, paste0(curname, ".markers.csv"))

  dms=tapply(markers$gene,markers$cluster,list)

  idy=1
  for(idy in c(1:length(dms))) {
    cluster = names(dms)[idy]
    allgenes=unlist(dms[idy])
    genes=unname(allgenes[1:min(gene_number, length(allgenes))])
    display_name=cnames[cluster,display_cluster_name]

    pdf(file=paste0(curname, ".", display_name, ".dot.pdf"), width=max(length(genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
    p<-DotPlot(subobj, group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
      xlab("genes") + ggtitle(display_name) + theme(plot.title = element_text(hjust = 0.5)) + xlab("Biomarker genes") + ylab(display_cluster_name)
    print(p)
    dev.off()
  }
}
