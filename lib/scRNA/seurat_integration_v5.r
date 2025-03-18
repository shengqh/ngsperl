rm(list=ls()) 
outFile='pbmc_rejection'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250226_T04_snRNA_hg38/seurat_rawdata/result/pbmc_rejection.rawobj.rds'
parFile2='/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250226_T04_snRNA_hg38/essential_genes/result/pbmc_rejection.txt'
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250226_T04_snRNA_hg38/seurat_sct2_rpca/result')

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
library(data.table)
library(SeuratData)
library(SeuratWrappers)
library(BiocParallel)

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

method=myoptions$integration_by_method_v5
reduction=myoptions$reduction

is_unix = .Platform$OS.type == "unix"
if(is_unix){
  ncores = as.numeric(myoptions$thread)
  bpparam = MulticoreParam(workers = ncores)
  register(bpparam)
}else{
  bpparam = SerialParam()
}

pca_dims<-1:as.numeric(myoptions$pca_dims)
cur_assay=ifelse(by_sctransform, "SCT", "RNA")

obj<-readRDS(parFile1)

# No matter scTransform or not, we need to normalize the object in order to get average expression later.
obj <- NormalizeData(obj, assay="RNA")

cat("preprocessing_rawobj ...\n")
finalList<-preprocessing_rawobj(
  rawobj=obj, 
  myoptions=myoptions, 
  prefix=detail_prefix, 
  filter_config_file=parFile3)

obj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

integrated_obj_file<-paste0(outFile, ".integrated.rds")
if(!file.exists(integrated_obj_file)){
  normalized_obj_file<-paste0(outFile, ".normalized.rds")
  if(!file.exists(normalized_obj_file)){
    if(has_batch_file){
      cat("Setting batch ...\n")
      poolmap = get_batch_samples(parSampleFile2, unique(obj$sample))
      obj@meta.data$batch <- unlist(poolmap[obj$sample])
    }else if(!("batch" %in% colnames(obj@meta.data))){
      obj@meta.data$batch <- obj$sample
    }

    #When using Seurat v5 assays, we can instead keep all the data in one object, but simply split the layers. 
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

    #integration would be done on the RNA assay
    DefaultAssay(obj) <- "RNA"

    if(by_sctransform){
      cat("SCTransform ...\n")
      obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)

      cat("JoinLayers of RNA assay ... \n")
      obj <- JoinLayers(obj, assay="RNA")

      DefaultAssay(obj) <- "SCT"
    }else{
      cat("NormalizeData/FindVariableFeatures ...\n")
      obj <- NormalizeData(obj)

      cat("FindVariableFeatures ... \n")
      obj <- FindVariableFeatures(obj)

      cat("ScaleData ... \n")
      obj <- ScaleData(obj, verbose = FALSE)
    }
    cat("Saving normalized object to file:", normalized_obj_file, "\n")
    saveRDS(obj, normalized_obj_file)
  }else{
    cat("Reading normalized object from file:", normalized_obj_file, "\n")
    obj<-readRDS(normalized_obj_file)
  }

  cat("RunPCA ... \n")
  obj <- RunPCA(object = obj, assay=cur_assay, verbose=FALSE)

  output_ElbowPlot(obj, detail_prefix, "pca")

  normalization_method = ifelse(by_sctransform, "SCT", "LogNormalize")

  if(method == "FastMNNIntegration"){
    cat("IntegrateLayers by FastMNNIntegration with thread", myoptions$thread, "... \n")
    obj <- IntegrateLayers(
      object = obj,
      method = FastMNNIntegration,
      orig.reduction = "pca",
      assay = cur_assay,
      new.reduction = reduction,
      verbose = T,
      #for fastMNN
      batch = obj@meta.data$batch, #batch is required for integration when only one object provided
      BPPARAM = bpparam
    )
  }else if (method == "CCAIntegration" | method == "RPCAIntegration"){
    obj <- IntegrateLayers(
      object = obj,
      method = method,
      orig.reduction = "pca",
      assay = cur_assay,
      new.reduction = reduction,
      normalization.method = normalization_method,
      verbose = T
    )
  }else{
    obj <- IntegrateLayers(
      object = obj,
      method = method,
      orig.reduction = "pca",
      assay = cur_assay,
      new.reduction = reduction,
      verbose = T
    )
  }

  #The FastMNNIntegration will create a new Seurat object with cur_array (might be SCT).
  #The default assay of new object will be set to RNA. It will cause problem when we use FindNeighbors and FindClusters.
  #So we need to set the default assay to cur_assay.
  obj[[reduction]]@assay.used = cur_assay

  cat("FindNeighbors ... \n")
  obj <- FindNeighbors( obj, 
                        dims = 1:30, 
                        assay=cur_assay, 
                        reduction = reduction)

  cat("RunUMAP ... \n")
  obj <- RunUMAP( obj, 
                  assay = cur_assay,
                  reduction = reduction, 
                  dims = 1:30)    

  if("ADT" %in% names(obj)){
    obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
  }
  cat("Saving integrated object to file:", integrated_obj_file, "\n")
  saveRDS(obj, integrated_obj_file)

  unlink(normalized_obj_file)
}else{
  cat("Reading integrated object from file:", integrated_obj_file, "\n")
  obj<-readRDS(integrated_obj_file)
}

DefaultAssay(obj)<-cur_assay
finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
cat("Saving final object to file:", finalListFile, "\n")
saveRDS(finalList, file=finalListFile)

rm(finalList)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj, detail_prefix, FALSE, myoptions$qc_genes)

unlink(integrated_obj_file)
