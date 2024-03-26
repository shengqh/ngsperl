rm(list=ls()) 
outFile='crs'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/h_turner_lab/shengq2/20221206_7114_8822_scRNA_hg38/seurat_rawdata/result/crs.rawobj.rds'
parFile2='/nobackup/h_turner_lab/shengq2/20221206_7114_8822_scRNA_hg38/essential_genes/result/crs.txt'
parFile3=''


setwd('/nobackup/h_turner_lab/shengq2/20221206_7114_8822_scRNA_hg38/seurat_sct_integration/result')

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

options(future.globals.maxSize= 10779361280)
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

cat("preprocessing_rawobj ...\n")
finalList<-preprocessing_rawobj(rawobj, myoptions, detail_prefix)
rawobj<-finalList$rawobj
finalList$rawobj<-NULL

if(has_batch_file){
  cat("Setting batch ...\n")
  poolmap = get_batch_samples(parSampleFile2, unique(rawobj$sample))
  rawobj$batch <- unlist(poolmap[rawobj$sample])
}else if(!("batch" %in% colnames(rawobj@meta.data))){
  rawobj$batch <- rawobj$sample
}

objs<-SplitObject(object = rawobj, split.by = "batch")
rm(rawobj)

if(length(objs) == 1){
  obj <- objs[[1]]
  if(by_sctransform){
    cat("SCTransform ...\n")
    obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
  }
}else{
  if(by_sctransform){
    # https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
    cat("SCTransform ...\n")
    # perform sctransform
    objs <- lapply(objs, FUN=SCTransform, method = "glmGamPoi")

    cat("SelectIntegrationFeatures ...\n")
    features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
    assay="SCT"
    all_genes <- Reduce(intersect, lapply(objs, rownames)) 

    cat("PrepSCTIntegration ...\n")
    objs <- PrepSCTIntegration(object.list = objs, anchor.features = features, verbose = FALSE)

    cat("FindIntegrationAnchors ...\n")
    obj.anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT", anchor.features = features, verbose = FALSE)

    cat("IntegrateData ...\n")
    obj <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", verbose = FALSE, features.to.integrate=all_genes)      
    DefaultAssay(obj) <- "integrated"  
  }else{
    # https://satijalab.org/seurat/articles/integration_introduction.html
    cat("NormalizeData/FindVariableFeatures ...\n")
    # normalize and identify variable features for each dataset independently
    objs<-lapply(objs, function(x){
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)  
      return(x)
    })  
    assay="RNA"
    all_genes <- Reduce(intersect, lapply(objs, rownames)) 

    cat("SelectIntegrationFeatures ...\n")
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = objs)

    cat("FindIntegrationAnchors ...\n")
    obj.anchors <- FindIntegrationAnchors(object.list = objs, dims = 1:20, anchor.features = features)

    # this command creates an 'integrated' data assay
    cat("IntegrateData ...\n")
    obj <- IntegrateData(anchorset = obj.anchors, dims = 1:20, features.to.integrate=all_genes)  

    # specify that we will perform downstream analysis on the corrected data note that the
    # original unmodified data still resides in the 'RNA' assay
    DefaultAssay(obj) <- "integrated"  
  }
}
rm(objs)

cat("ScaleData ... \n")
obj <- ScaleData(obj, verbose = FALSE)

cat("RunPCA ... \n")
obj <- RunPCA(object = obj, verbose=FALSE)

output_ElbowPlot(obj, detail_prefix, "pca")

cat("RunUMAP ... \n")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)    

if("ADT" %in% names(obj)){
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
}

finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
saveRDS(finalList, file=finalListFile)

rm(finalList)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj, detail_prefix, FALSE, myoptions$qc_genes)
