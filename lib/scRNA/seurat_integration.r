rm(list=ls()) 
outFile='PH_combine'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/cqs/shengq2/paula_hurley_projects/20221115_scRNA_7467_benign_hg38/seurat_rawdata/result/PH_combine.rawobj.rds'
parFile2='/scratch/cqs/shengq2/paula_hurley_projects/20221115_scRNA_7467_benign_hg38/essential_genes/result/PH_combine.txt'
parFile3=''


setwd('/scratch/cqs/shengq2/paula_hurley_projects/20221115_scRNA_7467_benign_hg38/seurat_sct_integration/result')

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

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

prefix<-outFile

has_batch_file<-file.exists(parSampleFile2)

pca_dims<-1:as.numeric(myoptions$pca_dims)

finalList<-list()
finalListFile<-paste0(prefix, ".preprocess.rds")

rawobj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(rawobj, myoptions, prefix)
rawobj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

poolmap = get_batch_samples(parSampleFile2, unique(rawobj$sample))

rawobj$batch=unlist(poolmap[rawobj$sample])
objs<-SplitObject(object = rawobj, split.by = "batch")
rm(rawobj)

#filter cells
finalList$filter<-list(nFeature_cutoff_min=myoptions$nFeature_cutoff_min,
                      nFeature_cutoff_max=myoptions$nFeature_cutoff_max,
                      mt_cutoff=myoptions$mt_cutoff,
                      nCount_cutoff=myoptions$nCount_cutoff,
                      nCount_sd_cutoff=myoptions$nCount_sd_cutoff)

if(length(objs) == 1){
  obj <- objs[[1]]
  if(by_sctransform){
    obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
  }
}else{
  if(by_sctransform){
    # https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
    cat("performing SCTransform ...")
    # perform sctransform
    objs <- lapply(objs, FUN=SCTransform, method = "glmGamPoi")
    features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
    assay="SCT"
    all_genes <- Reduce(intersect, lapply(objs, rownames)) 
    objs <- PrepSCTIntegration(object.list = objs, anchor.features = features, verbose = FALSE)
    obj.anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
    obj <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", verbose = FALSE, features.to.integrate=all_genes)      
    DefaultAssay(obj) <- "integrated"  
  }else{
    # https://satijalab.org/seurat/articles/integration_introduction.html
    cat("performing NormalizeData/FindVariableFeatures ...")

    # normalize and identify variable features for each dataset independently
    objs<-lapply(objs, function(x){
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)  
      return(x)
    })  
    assay="RNA"
    all_genes <- Reduce(intersect, lapply(objs, rownames)) 

    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = objs)

    obj.anchors <- FindIntegrationAnchors(object.list = objs, dims = 1:20, anchor.features = features)

    # this command creates an 'integrated' data assay
    obj <- IntegrateData(anchorset = obj.anchors, dims = 1:20, features.to.integrate=all_genes)  

    # specify that we will perform downstream analysis on the corrected data note that the
    # original unmodified data still resides in the 'RNA' assay
    DefaultAssay(obj) <- "integrated"  
  }
}

cat("ScaleData ... \n")
obj <- ScaleData(obj, verbose = FALSE)

cat("RunPCA ... \n")
obj <- RunPCA(object = obj, verbose=FALSE)

output_ElbowPlot(obj, outFile, "pca")

cat("RunUMAP ... \n")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)    

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj, outFile, FALSE, myoptions$qc_genes)

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
