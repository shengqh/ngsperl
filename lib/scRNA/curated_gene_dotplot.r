source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(ggpubr)
library(dendextend)
library(ggdendro)

finalList<-readRDS(parFile1)
obj<-finalList$obj

assay=ifelse("SCT" %in% names(obj@assays), "SCT", "RNA")

genes_df=read.table(parSampleFile1, sep="\t", stringsAsFactors=F)

cluster_df=read.table(parSampleFile2, sep="\t", stringsAsFactors=F)
cluster_request <- tapply(cluster_df$V1,cluster_df$V2,list)

params_def=read.table(parSampleFile3, sep="\t", stringsAsFactors=F)
params=tapply(params_def$V1,params_def$V2,list)

cell_df<-read_cell_cluster_file(parFile2)
cell_df<-read_cell_cluster_file(parFile2, sort_cluster_name = params$sort_cluster_name, display_cluster_name=params$display_cluster_name)

obj[["final_seurat_clusters"]]=cell_df[,params$display_cluster_name]

#assay=ifelse(params$by_sctransform=="1", "SCT", "RNA")
assaydata=GetAssayData(obj, assay=assay)
allgenes=rownames(assaydata)
rm(assaydata)
genes_df=subset(genes_df, genes_df$V1 %in% allgenes)

hasGroup = exists("parSampleFile4")
if(hasGroup){
  groups=read.table(parSampleFile4, sep="\t")
  gmap=split(groups$V2, groups$V1)
  obj$group=unlist(gmap[unlist(obj$orig.ident)])
  groupnames=unique(obj$group)
}

drawDotPlot<-function(object, geneset_name, genes, assay, split.by=NA){
  p<-DotPlot(object, assay = assay, group.by="final_seurat_clusters", features=genes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    theme(plot.title = element_text(hjust = 0.5)) + xlab(gsub("_", " ", geneset_name)) + ylab("")
  cat(geneset_name, "sortName \n")
  
  valid_genes=unique(p$data$features.plot)
  valid_genes=valid_genes[order(as.character(valid_genes))]
  width=max(length(valid_genes) * 0.3 + 2, 10)

  p$data$features.plot<-factor(p$data$features.plot, levels=valid_genes)
  pdata<-p$data
  
  cl=length(unlist(unique(object[["final_seurat_clusters"]])))
  height=max(8, cl * 0.3 + 1)
  
  pdf(file=paste0(geneset_name, ".dot.sortName.pdf"), width=width, height=height)
  print(p)
  dev.off()

  cat(geneset_name, "sortPct \n")

  adata<-acast(pdata, features.plot~id, value.var = "pct.exp")
  hc <- hclust(dist(adata))
  ordered_genes=rownames(adata)[order.hclust(hc)]
  p$data$features.plot<-factor(p$data$features.plot, ordered_genes)

  pdf(file=paste0(geneset_name, ".dot.sortPct.pdf"), width=width, height=height)
  print(p)
  dev.off()

  cat(geneset_name, "sortExp \n")
  
  minvalue=min(pdata$avg.exp.scaled[!is.na(pdata$avg.exp.scaled)])
  pdata$avg.exp.scaled[is.na(pdata$avg.exp.scaled)] = minvalue

  adata<-acast(pdata, features.plot~id, value.var = "avg.exp.scaled")
  hc <- hclust(dist(adata))
  ordered_genes=rownames(adata)[order.hclust(hc)]
  p$data$features.plot<-factor(p$data$features.plot, ordered_genes)
  
  pdf(file=paste0(geneset_name, ".dot.sortExp.pdf"), width=width, height=height)
  print(p)
  dev.off()
  
  if(!is.na(split.by)){
    cat(geneset_name, "~", split.by, "sortName \n")
    #split.by="group"

    groupnames=unique(object[[split.by, drop=TRUE]])
    cols=rainbow(length(groupnames))
    names(cols)=groupnames
    
    p<-DotPlot(object, assay = assay, group.by="final_seurat_clusters", features=genes, split.by=split.by, dot.scale = 8) + RotatedAxis() +
      theme(plot.title = element_text(hjust = 0.5)) + xlab(gsub("_", " ", geneset_name)) + ylab("")
    
    valid_genes=unique(p$data$features.plot)
    valid_genes=valid_genes[order(as.character(valid_genes))]
    p$data$features.plot<-factor(p$data$features.plot, levels=valid_genes)
    
    data.plot=p$data
    data.plot$split.by=unlist(lapply(data.plot$id, function(x){
      res=x
      while(!(res %in% groupnames)){
        res=gsub(".+?_", "", x)
      }
      return(res)
    }))
    
    data.plot$split.color=cols[data.plot$split.by]
    
    data.plot$colors=unlist(apply(data.plot, 1, function(x){
      color=x[["split.color"]]
      value=as.numeric(x[["avg.exp.scaled"]])
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }))
    
    p$data$colors=data.plot$colors

    cl=length(unlist(unique(object[["final_seurat_clusters"]]))) * length(groupnames)
    
    pdf(file=paste0(geneset_name, ".dot.", split.by, ".sortName.pdf"), width=width, height=max(8, cl * 0.3 + 1))
    print(p)
    dev.off()
  }
}

idx=1
for(idx in c(1:length(cluster_request))) {
  geneset_name=names(cluster_request)[idx]
  clusternames = as.character(cluster_request[[idx]])

  if(clusternames=='all'){
    subobj=obj
  }else{
    cells=rownames(cell_df)[cell_df[,params$cluster_name] %in% clusternames]
    subobj=subset(obj, cells=cells)
  }

  cur_genes_df=subset(genes_df, genes_df$V2 == geneset_name)
  genes=unique(cur_genes_df$V1)

  if(hasGroup){
    drawDotPlot(subobj, geneset_name, genes, assay, split.by="group")
  }else{
    drawDotPlot(subobj, geneset_name, genes, assay)
  }
}
