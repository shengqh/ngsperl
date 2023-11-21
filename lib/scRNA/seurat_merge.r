rm(list=ls()) 
outFile='P5798'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/seurat_rawdata/result/P5798.rawobj.rds'
parFile2='/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/essential_genes/result/P5798.txt'
parFile3='/data/cqs/shengq2/program/collaborations/jonathan_brown/20231114_scRNA_5798_human_liver_redo/filter_config.csv'


setwd('/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/seurat_sct2_merge/result')

### Parameter setting end ###

source("scRNA_func.r")
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(cowplot)
library(scales)
library(stringr)
library(glmGamPoi)
require(data.table)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
npcs<-as.numeric(myoptions$pca_dims)
pca_dims<-1:npcs

by_sctransform<-is_one(myoptions$by_sctransform)
use_sctransform_v2<-is_one(myoptions$use_sctransform_v2)
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)
thread<-to_numeric(myoptions$thread, 1)
is_preprocessed=is_one(myoptions$is_preprocessed)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile2, sep="\t" ,header=F)$V1

prefix<-outFile

species=myoptions$species

finalListFile<-paste0(prefix, ".final.rds")

obj<-readRDS(parFile1)

if(!is_preprocessed){
  finalList<-preprocessing_rawobj(obj, myoptions, prefix, parFile3)
  obj<-finalList$rawobj
  finalList<-finalList[names(finalList) != "rawobj"]

  DefaultAssay(obj)<-"RNA"
  if(by_sctransform){
    obj<-do_sctransform(obj, vars.to.regress=vars.to.regress, return.only.var.genes=FALSE, mc.cores=thread, use_sctransform_v2=use_sctransform_v2)
  }
}else{
  finalList=list()
  if(by_sctransform){
    #https://github.com/satijalab/seurat/issues/2814
    VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data)
  }
}
assay=ifelse(by_sctransform, "SCT", "RNA")

#no matter if we will use sctransform, we need normalized RNA assay for visualization and cell type annotation
#data slot for featureplot, dotplot, cell type annotation and scale.data slot for heatmap
obj<-do_normalization(obj, selection.method="vst", nfeatures=3000, vars.to.regress=vars.to.regress, scale.all=FALSE, essential_genes=essential_genes)

DefaultAssay(obj)<-assay

cat("RunPCA ... \n")
obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)

# #https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# pcs <- find_number_of_reduction(obj, reduction="pca")
# final_pcs=max(pcs, as.numeric(myoptions$pca_dims))
# writeLines(paste0("pcs\t", final_pcs), con=paste0(outFile, ".pca.txt"))
# cat(paste0("recommend pcs=", final_pcs))

output_ElbowPlot(obj, outFile, "pca")

cat("RunUMAP ... \n")
obj <- RunUMAP(object = obj, dims=pca_dims, verbose = FALSE)

if(by_sctransform){
  #clear SCT data for small object
  #https://github.com/satijalab/seurat/issues/2587
  obj[["SCT"]]@misc <- NULL    
}    

finalList$obj<-obj

cat("saveRDS ... \n")
saveRDS(finalList, file=finalListFile)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj, outFile, FALSE, myoptions$qc_genes)
