rm(list=ls()) 
outFile='GPA_NML'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20230108_9112_3885_JH_scRNA/seurat_sct_harmony/result/GPA_NML.final.rds'
parFile2='/data/h_gelbard_lab/projects/20230108_9112_3885_JH_scRNA/seurat_sct_harmony_dr0.5_nrh_01_call/result/GPA_NML.scDynamic.meta.rds'
parFile3=''


setwd('/data/h_gelbard_lab/projects/20230108_9112_3885_JH_scRNA/seurat_sct_harmony_dr0.5_nrh_01_call_bubblemap/result')

### Parameter setting end ###

source("scRNA_func.r")
library("Seurat")
library("readxl")
library(ggplot2)

bubble_files=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
cluster_name = ifelse(is.null(myoptions$cluster_name), "seurat_clusters", myoptions$cluster_name)
celltype_name = ifelse(is.null(myoptions$celltype_name), "cell_type", myoptions$celltype_name)

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2)
}

bnames=unique(bubble_files$V3)
bn=bnames[1]
for(bn in bnames){
  cat(bn, "\n")
  cluster_pattern=bubble_files$V1[bubble_files$V2 == "cluster_pattern" & bubble_files$V3 == bn][1]
  bubblemap_file=bubble_files$V1[bubble_files$V2 == "file" & bubble_files$V3 == bn][1]
  ignore.case=is_one(bubble_files$V1[bubble_files$V2 == "ignore.case" & bubble_files$V3 == bn][1])
  rotate.title= is_one(bubble_files$V1[bubble_files$V2 == "rotate.title" & bubble_files$V3 == bn][1])
  width=as.numeric(bubble_files$V1[bubble_files$V2 == "width" & bubble_files$V3 == bn][1])
  height=as.numeric(bubble_files$V1[bubble_files$V2 == "height" & bubble_files$V3 == bn][1])
  
  if(cluster_pattern != "" & cluster_pattern != "*"){
    cells = colnames(obj)[grepl(cluster_pattern, obj$cell_type, ignore.case = ignore.case)]
    subobj=subset(obj, cells=cells)
  }else{
    subobj=obj
  }

  g=get_bubble_plot(subobj, cluster_name, celltype_name, bubblemap_file, assay="RNA", orderby_cluster=TRUE, rotate.title=rotate.title) + scale_color_gradient2(low="blue", mid="yellow", high="red")
  subobj$dump_cluster = paste0(unlist(subobj[[cluster_name]]), ":", unlist(subobj[[celltype_name]]))
  png(paste0(outFile, ".", bn, ".bubblemap.png"), width=width, height=get_dot_height(subobj, "dump_cluster"),res=300)
  print(g)
  dev.off()

  g=get_bubble_plot(subobj, cluster_name, "orig.ident", bubblemap_file, assay="RNA", orderby_cluster=TRUE, rotate.title=rotate.title) + scale_color_gradient2(low="blue", mid="yellow", high="red")
  subobj$dump_cluster = paste0(unlist(subobj[[cluster_name]]), ":", unlist(subobj[["orig.ident"]]))
  png(paste0(outFile, ".", bn, ".bubblemap.ct_ident.png"), width=width, height=get_dot_height(subobj, "dump_cluster"),res=300)
  print(g)
  dev.off()

  g=get_bubble_plot(subobj, "orig.ident", cluster_name, bubblemap_file, assay="RNA", orderby_cluster=TRUE, rotate.title=rotate.title) + scale_color_gradient2(low="blue", mid="yellow", high="red")
  subobj$dump_cluster = paste0(unlist(subobj[["orig.ident"]]), ":", unlist(subobj[[cluster_name]]))
  png(paste0(outFile, ".", bn, ".bubblemap.ident_ct.png"), width=width, height=get_dot_height(subobj, "dump_cluster"),res=300)
  print(g)
  dev.off()
}
