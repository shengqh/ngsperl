rm(list=ls()) 
outFile='PEO4'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/h_vangard_1/shengq2/guoyan/20231004_pipseq_10x_comparison/20231005_scRNA/seurat_rawdata/result/PEO4.rawobj.rds'
parFile2='/nobackup/h_vangard_1/shengq2/guoyan/20231004_pipseq_10x_comparison/20231005_scRNA/essential_genes/result/PEO4.txt'
parFile3=''


setwd('/nobackup/h_vangard_1/shengq2/guoyan/20231004_pipseq_10x_comparison/20231005_scRNA/seurat_sct2_harmony/result')

### Parameter setting end ###

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(cowplot)
library(scales)
library(stringr)
library(harmony)
library(patchwork)
require(data.table)

source("scRNA_func.r")

options(future.globals.maxSize=1024^3*100) #100G
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

detail_folder=paste0(myoptions$task_name, gsub(".html","",myoptions$rmd_ext))
dir.create(detail_folder, showWarnings = FALSE)
detail_prefix=file.path(detail_folder, myoptions$task_name)

has_batch_file<-file.exists(parSampleFile2)
npcs<-as.numeric(myoptions$pca_dims)
pca_dims<-1:npcs

obj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(obj, myoptions, detail_prefix)
obj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

essential_genes=read.table(parFile2, sep="\t" ,header=F)$V1

by_sctransform<-is_one(myoptions$by_sctransform)
use_sctransform_v2<-is_one(myoptions$use_sctransform_v2)

regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}


ignore_variable_genes=c()
if("ignore_variable_gene_file" %in% names(myoptions)){
  ignore_variable_gene_file = myoptions$ignore_variable_gene_file
  if(ignore_variable_gene_file != ""){
    ignore_variable_genes=readLines(ignore_variable_gene_file)
  }
}

obj<-do_harmony(
  obj = obj, 
  by_sctransform = by_sctransform, 
  vars.to. = vars.to.regress, 
  has_batch_file = has_batch_file, 
  batch_file = parSampleFile2, 
  pca_dims = pca_dims,
  essential_genes=essential_genes, 
  mc.cores=8, 
  use_sctransform_v2=use_sctransform_v2,
  ignore_variable_genes=ignore_variable_genes)

reduction="harmony"

for (reduct in c("pca", "harmony")){
  output_ElbowPlot(obj, detail_prefix, reduct)
}

cat("RunUMAP ... ")
obj <- RunUMAP(object = obj, reduction=reduction, dims=pca_dims, verbose = FALSE)

if("ADT" %in% names(obj)){
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
}

finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
saveRDS(finalList, file=finalListFile)

output_integration_dimplot(obj, detail_prefix, FALSE, myoptions$qc_genes)
