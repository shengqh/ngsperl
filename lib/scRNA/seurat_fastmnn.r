rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_rawdata/result/combined.rawobj.rds'
parFile2='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/essential_genes/result/combined.txt'
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_sct_fastmnn/result')

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
library(harmony)
library(patchwork)
require(data.table)
library(SeuratWrappers)
library(SeuratData)

options(future.globals.maxSize=1024^3*100) #100G
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

detail_folder=paste0(myoptions$task_name, gsub(".html","",myoptions$rmd_ext))
dir.create(detail_folder, showWarnings = FALSE)
detail_prefix=file.path(detail_folder, myoptions$task_name)

has_batch_file<-file.exists(parSampleFile2)

pca_dims<-1:as.numeric(myoptions$pca_dims)

rawobj<-readRDS(parFile1)

cat("preprocessing_obj ...\n")
finalList<-preprocessing_rawobj(
    rawobj=rawobj, 
    myoptions=myoptions, 
    prefix=detail_prefix, 
    filter_config_file=parFile3)
rawobj<-finalList$rawobj
finalList$rawobj<-NULL

if(has_batch_file){
  cat("Setting batch ...\n")
  poolmap = get_batch_samples(parSampleFile2, unique(rawobj$sample))
  rawobj$batch <- unlist(poolmap[rawobj$sample])
}else if(!("batch" %in% colnames(rawobj@meta.data))){
  rawobj$batch <- rawobj$sample
}

##### Run FastMNN #####
cat("NormalizeData ... \n")
rawobj <- NormalizeData(rawobj)

cat("FindVariableFeatures ... \n")
rawobj <- FindVariableFeatures(rawobj, nfeatures = 3000)

# cat("RunPCA ... \n")
# rawobj <- RunPCA(object = rawobj, verbose=FALSE)

cat("SplitObject ... \n")
objs<-SplitObject(object = rawobj, split.by = "batch")
rm(rawobj)

cat("RunFastMNN ... \n")
obj <- RunFastMNN(object.list = objs, features = 3000)
rm(objs)

cat("FindVariableFeatures ... \n")
obj <- FindVariableFeatures(obj, nfeatures = 3000)

cat("ScaleData ... \n")
obj <- ScaleData(obj, verbose = FALSE)

cat("RunPCA ... \n")
obj <- RunPCA(object = obj, verbose=FALSE)

output_ElbowPlot(obj, detail_prefix, "pca")

cat("RunUMAP ... \n")
obj <- RunUMAP(obj, reduction = "mnn", dims = 1:30)

cat("FindNeighbors ... \n")
obj <- FindNeighbors(obj, reduction = "mnn", dims = 1:30)

if("ADT" %in% names(obj)){
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
}

finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
saveRDS(finalList, file=finalListFile)

#finalList<-readRDS(finalListFile)
#obj=finalList$obj

rm(finalList)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj=obj, 
  outFile=detail_prefix, 
  has_batch_file=FALSE, 
  qc_genes=myoptions$qc_genes)

