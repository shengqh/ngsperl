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
library(future)

options(future.globals.maxSize= 10779361280)
future::plan(strategy = "multisession")

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

thread=as.numeric(myoptions$thread)
is_unix = .Platform$OS.type == "unix"
if(is_unix){
  bpparam = MulticoreParam(workers = thread)
  register(bpparam)
}else{
  bpparam = SerialParam()
}

pca_dims<-1:as.numeric(myoptions$pca_dims)
cur_assay=ifelse(by_sctransform, "SCT", "RNA")

ignore_variable_genes=c()
if("ignore_variable_gene_file" %in% names(myoptions)){
  ignore_variable_gene_file = myoptions$ignore_variable_gene_file
  if(ignore_variable_gene_file != ""){
    ignore_variable_genes=readLines(ignore_variable_gene_file)
  }
}

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

    #integration would be done on the RNA assay
    DefaultAssay(obj) <- "RNA"

    # In order to perform integration, we need to split the object by batch, no matter through SCTransform or not.
    # When using Seurat v5 assays, we can instead keep all the data in one object, but simply split the layers. 
    cat("Split RNA assay ...\n")
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)

    if(by_sctransform){
      cat("SCTransform ...\n")
      obj <- SCTransform(obj, method = "glmGamPoi", verbose = FALSE)
      DefaultAssay(obj) <- "SCT"

      cat("JoinLayers of RNA assay ... \n")
      obj <- JoinLayers(obj, assay="RNA")
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

  if(length(ignore_variable_genes) > 0){
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), ignore_variable_genes)
  }

  obj = do_PCA_Integration( obj, 
                            cur_assay, 
                            by_sctransform, 
                            method=method, 
                            new.reduction=reduction, 
                            orig.reduction="pca",
                            thread=thread,
                            detail_prefix=detail_prefix)

  if(!by_sctransform){
    cat("JoinLayers of", cur_assay, "assay ... \n")
    obj <- JoinLayers(obj, assay=cur_assay)
  }
  
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
