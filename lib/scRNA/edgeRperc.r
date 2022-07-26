rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.final.rds'
parFile2='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample/result/AK6383.edgeR.files.csv'
parFile3='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.meta.rds'


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_edgeR_inCluster_bySample_vis/result')

### Parameter setting end ###

source("scRNA_func.r")
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(R.utils)
library(reshape2)
library(patchwork)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
cluster_name=myoptions$cluster_name

if(!exists("obj")){
  obj<-read_object(parFile1, parFile3, cluster_name)
}

g<-DimPlot(obj, group.by = cluster_name, label=T) + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5))

edgeRres<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
edgeRfolder<-dirname(parFile2)
rownames(edgeRres)<-edgeRres$prefix

df<-data.frame(c1=obj$seurat_clusters, c2=obj[[cluster_name]])
df<-unique(df)
df<-df[order(df$c1),]
obj[[cluster_name]]<-factor(unlist(obj[[cluster_name]]), levels=unique(df[,cluster_name]))

result<-NULL
prefix<-rownames(edgeRres)[4]
for (prefix in rownames(edgeRres)){
  cat("Processing ", prefix, "\n")
  comparison<-edgeRres[prefix, "comparison"]
  sigGenenameFile<-paste0(edgeRfolder, "/", edgeRres[prefix, "sigGenenameFile"])
  cellType<-edgeRres[prefix, "cellType"]
  deFile=gsub(".sig_genename.txt", ".csv", sigGenenameFile)
  totalGene=length(readLines(deFile))-1
  sigGene=length(readLines(sigGenenameFile))
  curDF<-data.frame("prefix"=prefix, "sigGene"=sigGene, "totalGene"=totalGene, "cluster"=cellType, "comparison"=comparison)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

result$sigRate<-result$sigGene * 100.0 / result$totalGene

comp=unique(result$comparison)[1]
for (comp in unique(result$comparison)) {
  compRes = result[result$comparison == comp,]
  rateMap=unlist(split(compRes$sigRate, compRes$cluster))
  
  obj$sigRate=rateMap[as.character(unlist(obj[[cluster_name]]))]
  gf<-MyFeaturePlot(obj, feature="sigRate", cols=c("lightgrey", "red")) + ggtitle(comp) + theme(plot.title = element_text(hjust=0.5))
  g<-g+gf
}

n<-length(unique(result$comparison))+1
ncol<-ceiling(sqrt(n))

g<-g+plot_layout(ncol=ncol)
png(paste0(outFile, ".sigGenePerc.png"), width=3300, height=3000, res=300)
print(g)
dev.off()
