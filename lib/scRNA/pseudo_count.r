rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_05_pseudo_count/result')

### Parameter setting end ###

source("scRNA_func.r")
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(testit)
library(tools)
library(Matrix.utils)

sumcount<-function(ct_count, groupings){
  rescount<-t(aggregate.Matrix(t(ct_count), groupings=groupings, fun="sum"))
  psum<-rowSums(rescount)
  rescount<-rescount[psum>0,]
  return(rescount)
}

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
cluster_name=myoptions$group_by

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2, cluster_name)
}

clusterDf<-obj@meta.data

cts = as.character(unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name]))
prefixList<-gsub('[ /:_]+', '_', cts)

res_files=c()

idx<-1
for (idx in c(1:length(cts))){
  ct = cts[idx]
  cat(ct, "\n")
  
  prefix = prefixList[idx]
  
  clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
  de_obj<-subset(obj, cells=rownames(clusterCt))

  ct_count<-de_obj@assays$RNA@counts
  groupings<-unlist(de_obj$orig.ident)
  p_count<-sumcount(ct_count, groupings)

  p_file=paste0(prefix, ".pseudo_count.csv")
  
  write.csv(p_count, p_file)
  
  res_files<-c(res_files, file_path_as_absolute(p_file))
}    

res_df<-data.frame("cluster"=cts, "prefix"=prefixList, "pusedo_file"=p_file)
write.csv(res_df, paste0(outFile, ".pusedo_count.list.csv"), row.names=F)


