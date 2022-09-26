#rm(list=ls()) 
outFile='PL_7114_human'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='C:/projects/nobackup/h_turner_lab/shengq2/20220805_7114_scRNA_human/seurat_sct_harmony_multires_03_choose/result/PL_7114_human.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/nobackup/h_turner_lab/shengq2/20220805_7114_scRNA_human/seurat_sct_harmony_multires_03_choose_bubblemap/result')

### Parameter setting end ###

source("scRNA_func.r")
library("Seurat")
library("readxl")
library(ggplot2)

bubble_files=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-unlist(split(options_table$V1, options_table$V2))

if(!exists('obj')){
  obj<-read_object(parFile1)
}

bnames=unique(bubble_files$V3)
bn=bnames[3]
for(bn in bnames){
  cat(bn, "\n")
  cluster_pattern=bubble_files$V1[bubble_files$V2 == "cluster_pattern" & bubble_files$V3 == bn][1]
  bubblemap_file=bubble_files$V1[bubble_files$V2 == "file" & bubble_files$V3 == bn][1]
  ignore.case=is_one(bubble_files$V1[bubble_files$V2 == "ignore.case" & bubble_files$V3 == bn][1])
  rotate.title= is_one(bubble_files$V1[bubble_files$V2 == "rotate.title" & bubble_files$V3 == bn][1])
  
  if(cluster_pattern == "" | cluster_pattern == "*"){
    subobj=obj
  }else{
    cells = colnames(obj)[grepl(cluster_pattern, obj$cell_type, ignore.case = ignore.case)]
    subobj=subset(obj, cells=cells)
  }
  g=get_bubble_plot(subobj, "seurat_clusters", "cell_type", bubblemap_file, assay="RNA", orderby_cluster=TRUE, rotate.title=rotate.title)
  
  png(paste0(outFile, ".", bn, ".bubblemap.png"), width=5000, height=2000,res=300)
  print(g)
  dev.off()
}

