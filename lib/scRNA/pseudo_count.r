rm(list=ls()) 
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/cellbender_nd_seurat_fastmnn_dr0.5_3_choose/result/Aorta_Progeria.final.rds'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/cellbender_nd_seurat_fastmnn_dr0.5_3_choose_pesudo_count/result')

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
cluster_name=ifelse(is.null(myoptions$group_by), "seurat_cell_type", myoptions$group_by)
sample_name=ifelse(is.null(myoptions$sample_name), "orig.ident", myoptions$sample_name)

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2, cluster_name)
}

clusterDf<-obj@meta.data

cts = as.character(unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name]))
prefixList<-celltype_to_filename(cts)

res_files=c()

cat("Get bulk pesudo count\n")
bulk_count=get_group_count(obj, sample_name) 
bulk_count_file=paste0(outFile, ".pseudo_count.bulk.csv")
write.csv(bulk_count, file = bulk_count_file, quote = F)

res_files=c(res_files, file_path_as_absolute(bulk_count_file))

idx<-1
for (idx in c(1:length(cts))){
  ct = cts[idx]
  cat("Get", ct, "pesudo count\n")
  
  prefix = prefixList[idx]
  
  clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
  de_obj<-subset(obj, cells=rownames(clusterCt))

  bulk_count=get_group_count(de_obj, sample_name) 
  bulk_count_file=paste0(outFile, ".pseudo_count", ".", prefix, ".csv")
  write.csv(bulk_count, file = bulk_count_file, quote = F)

  res_files=c(res_files, file_path_as_absolute(bulk_count_file))
}    

res_df<-data.frame("cluster"=c("bulk", cts), "prefix"=c("bulk", prefixList), "pseudo_file"=res_files)
write.csv(res_df, paste0(outFile, ".pseudo_count.list.csv"), row.names=F)
