rm(list=ls()) 
outFile='P9061'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/vickers_lab/projects/20251128_9061_scRNA_mm10_cellbender_redo/cellbender_scDblFinder_seurat_rawdata/result/P9061.rawobj.rds'
parFile2='/nobackup/vickers_lab/projects/20251128_9061_scRNA_mm10_cellbender_redo/essential_genes/result/P9061.txt'
parFile3=''


setwd('/nobackup/vickers_lab/projects/20251128_9061_scRNA_mm10_cellbender_redo/cellbender_scDblFinder_seurat_fastmnn/result')

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

options(future.globals.maxSize=1024^3*100) #100G

#Sometimes, using multisession for parallel processing would cause unexpected error. comment it out at 202505223.
#future::plan(strategy = "multisession")

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

cat("preprocessing_rawobj ...\n")
finalList<-preprocessing_rawobj(
  rawobj=obj, 
  myoptions=myoptions, 
  prefix=detail_prefix, 
  filter_config_file=parFile3)
cat("preprocessing_rawobj done.\n")

obj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

if(has_batch_file){
  cat("Setting batch ...\n")
  poolmap = get_batch_samples(parSampleFile2, unique(obj$sample))
  obj@meta.data$batch <- unlist(poolmap[obj$sample])
}else if(!("batch" %in% colnames(obj@meta.data))){
  obj@meta.data$batch <- obj$sample
}

obj <- do_integration_v5(
  outFile=outFile,
  obj=obj,
  by_sctransform=by_sctransform,
  cur_assay=cur_assay,
  method=method,
  reduction=reduction,
  thread=thread,
  detail_prefix=detail_prefix,
  ignore_variable_genes=ignore_variable_genes
)

DefaultAssay(obj)<-cur_assay
finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
cat("Saving final object to file:", finalListFile, "\n")
saveRDS(finalList, file=finalListFile)

rm(finalList)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(
  obj=obj, 
  outFile=detail_prefix, 
  has_batch_file=FALSE, 
  qc_genes=myoptions$qc_genes, 
  assay="RNA")

unlink(integrated_obj_file)

