rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20240325_scRNA_iSGS_cell_atlas/06_scRNA_nosct/seurat_rawdata/result/iSGS_cell_atlas.rawobj.rds'
parFile2='/data/h_gelbard_lab/projects/20240325_scRNA_iSGS_cell_atlas/06_scRNA_nosct/essential_genes/result/iSGS_cell_atlas.txt'
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240325_scRNA_iSGS_cell_atlas/06_scRNA_nosct/seurat_merge/result')

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

detail_folder=paste0(myoptions$task_name, gsub(".html","",myoptions$rmd_ext))
dir.create(detail_folder, showWarnings = FALSE)
detail_prefix=file.path(detail_folder, myoptions$task_name)

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

ignore_variable_genes=c()
if("ignore_variable_gene_file" %in% names(myoptions)){
  ignore_variable_gene_file = myoptions$ignore_variable_gene_file
  if(ignore_variable_gene_file != ""){
    ignore_variable_genes=readLines(ignore_variable_gene_file)
  }
}

finalListFile<-paste0(prefix, ".final.rds")

obj<-readRDS(parFile1)

if(!is_preprocessed){
  finalList<-preprocessing_rawobj(
      rawobj=obj, 
      myoptions=myoptions, 
      prefix=detail_prefix, 
      filter_config_file=parFile3)

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
obj<-do_normalization(obj=obj, 
                      selection.method="vst", 
                      nfeatures=2000, 
                      vars.to.regress=vars.to.regress, 
                      scale.all=FALSE, 
                      essential_genes=essential_genes,
                      ignore_variable_genes=ignore_variable_genes)

DefaultAssay(obj)<-assay

cat("RunPCA ... \n")
obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)

# #https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# pcs <- find_number_of_reduction(obj, reduction="pca")
# final_pcs=max(pcs, as.numeric(myoptions$pca_dims))
# writeLines(paste0("pcs\t", final_pcs), con=paste0(outFile, ".pca.txt"))
# cat(paste0("recommend pcs=", final_pcs))

output_ElbowPlot(obj, detail_prefix, "pca")

cat("RunUMAP ... \n")
obj <- RunUMAP(object = obj, dims=pca_dims, verbose = FALSE)

if(by_sctransform){
  #clear SCT data for small object
  #https://github.com/satijalab/seurat/issues/2587
  obj[["SCT"]]@misc <- NULL    
}    

if("ADT" %in% names(obj)){
  cat("NormalizeData of ADT ... \n")
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
}

finalList$obj<-obj

cat("saveRDS ... \n")
saveRDS(finalList, file=finalListFile)
#finalList<-readRDS(finalListFile)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj, detail_prefix, FALSE, myoptions$qc_genes)

