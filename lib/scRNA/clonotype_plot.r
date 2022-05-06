library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)
obj<-finalList$obj

assay=ifelse("SCT" %in% names(obj@assays), "SCT", "RNA")

cell_df<-read_cell_cluster_file(parFile2)

cluster_df=read.table(parSampleFile1)
cluster_request <- tapply(cluster_df$V1,cluster_df$V2,list)

params_def=read.table(parSampleFile2)
params <- setNames(as.character(params_def$V1), params_def$V2)
gene_number=as.numeric(params['gene_number'])
cluster_name=params['cluster_name']
display_cluster_name=params['display_cluster_name']

if(display_cluster_name == "seurat_clusters"){
  display_cluster_name = "display_seurat_clusters"
}
obj[["final_seurat_clusters"]]=cell_df[,display_cluster_name]

cnames=unique(cell_df[,c(cluster_name, display_cluster_name)])
rownames(cnames)=cnames[,cluster_name]
cnames[,"filename"] = gsub('[ :]+', '_', cnames[,display_cluster_name])

idx=1
for(idx in c(1:length(cluster_request))) {
  curname=names(cluster_request)[idx]
  clusternames = as.character(cluster_request[[idx]])

  cells=rownames(cell_df)[cell_df[,cluster_name] %in% clusternames]
  subobj=subset(obj, cells=cells)

  cname_df=data.frame("cluster_name"=subobj[[cluster_name]], "display_cluster_name"=subobj[["final_seurat_clusters"]])

  markers=FindAllMarkers(subobj, assay="RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  # Try to find pairwise biomarker
  # allmarkers=NULL
  # for(c1 in clusternames){
  #   c1markers=subset(markers, markers$cluster==c1)
  #   c1markers$sig_count=0
  #   for(c2 in clusternames){
  #     if(c1 == c2){
  #       next
  #     }
  #     c1c2_markers = FindMarkers(subobj, ident.1 = c1, ident.2=c2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #     sig_c1c2 = ifelse(rownames(c1markers) %in% rownames(c1c2_markers), 1, 0)
  #     c1markers$sig_count = c1markers$sig_count + sig_c1c2
  #   }
  #   allmarkers=rbind(allmarkers, c1markers)
  # }
  # allmarkers=subset(allmarkers, allmarkers$p_val_adj < 0.01)
  # allmarkers=allmarkers[order(allmarkers$cluster, -allmarkers$sig_count, allmarkers$p_val_adj),]
  # markers=allmarkers

  write.csv(markers, paste0(curname, ".markers.csv"))

  dms=tapply(markers$gene,markers$cluster,list)

  all_genes=c()
  idy=1
  for(idy in c(1:length(dms))) {
    cluster = names(dms)[idy]
    allgenes=unlist(dms[idy])
    genes=unname(allgenes[1:min(gene_number, length(allgenes))])
    display_name=as.character(cnames[cluster,display_cluster_name])
    file_name=cnames[cluster,"filename"]
    
    all_genes = c(all_genes, genes)

    pdf(file=paste0(curname, ".", file_name, ".dot.pdf"), width=max(length(genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
    p<-DotPlot(subobj, assay="RNA", group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
      xlab("") + ylab("")
    print(p)
    dev.off()
  }
  
  pdf(file=paste0(curname, ".dot.pdf"), width=max(length(all_genes) * 0.4, 10), height=max(6, min(10, length(clusternames))))
  p<-DotPlot(subobj, assay="RNA", group.by="final_seurat_clusters", features=all_genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab("") + ylab("") + ggtitle("Seurat Marker Genes") + theme(plot.title = element_text(hjust = 0.5))
  print(p)
  dev.off()
}
