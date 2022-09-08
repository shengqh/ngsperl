rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.final.rds'
parFile2='/home/shengq2/program/projects/alexander_gelbard/20220103_cell_markers.txt'
parFile3='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.meta.rds'


setwd('/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_genes_marker/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(patchwork)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

celltype_name=myoptions$celltype_name
cluster_name=myoptions$cluster_name
samples=myoptions$samples

obj<-read_object(parFile1)
geneFile<-parFile2

DefaultAssay(obj)<-"RNA"

celltypes<-read.table(geneFile, sep="\t", header=T, stringsAsFactors = F)
celltypes$Gene<-gsub("\\s.+", "", celltypes$Gene)
celltypes$Gene<-toupper(celltypes$Gene)

missGenes<-celltypes[!(celltypes$Gene %in% rownames(obj)),]
write.csv(missGenes, "miss_gene.csv", row.names=F)
celltypes<-celltypes[celltypes$Gene %in% rownames(obj),]

allgenes<-rownames(obj)
allgenes<-allgenes[order(allgenes)]
writeLines(allgenes, con="all_genes.txt")

if(grepl('.rds$', parFile3)){
  clusterDf<-readRDS(parFile3)
}else{
  clusterDf<-read.csv(parFile3, stringsAsFactors = F, row.names=1)
}

cts=cluster_to_cell_type(clusterDf)

cts<-sort_cell_type(cts, celltype_name)

clusterDf<-clusterDf[colnames(obj),]
clusterDf$sample<-obj$orig.ident

obj$final_seurat_clusters<-factor(clusterDf[,cluster_name], levels=cts[,cluster_name])

if (!is.null(samples) && samples !=""){
  sample_idents = unlist(strsplit(samples, ";"))
  sample_cells<-rownames(clusterDf)[clusterDf$sample %in% sample_idents]
  obj=subset(obj, cells=sample_cells)
}

ct<-unique(celltypes$Celltype)[1]
for(ct in unique(celltypes$Celltype)){
  cat(ct, "\n")
  subcelltypes<-celltypes[celltypes$Celltype==ct,]
  expressedGenes<-toupper(unique(subcelltypes$Gene[subcelltypes$Status=="expressed"]))
  absentGenes<-toupper(unique(subcelltypes$Gene[subcelltypes$Status=="absent"]))
  
  width=max(5000, nrow(subcelltypes) * 70)
  png(filename=paste0(outFile, ".", ct, ".png"), width=width, height=2500, res=300)
  p<-DotPlot(obj, assay = "RNA", group.by="final_seurat_clusters", features=expressedGenes, cols = c("lightgrey", "red"), dot.scale = 8) + RotatedAxis() +
    xlab(paste0(ct, " expressed genes"))
  if(length(absentGenes) > 0){
    e<-DotPlot(obj, assay = "RNA", group.by="final_seurat_clusters", features=absentGenes, cols = c("lightgrey", "blue"), dot.scale = 8) + RotatedAxis()+
      xlab(paste0(ct, " absent genes"))
    p2<-p+e
    print(p2)
  }else{
    print(p)
  }
  dev.off()
}
