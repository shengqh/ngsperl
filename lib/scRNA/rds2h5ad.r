rm(list=ls()) 
sample_name='MTC16'
outFile='MTC16'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/wenda_ye_projects/20261120_MTC_UAB_TDMI/CellTypist_1_h5ad/result/MTC16')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(SingleCellExperiment)
library(rhdf5)
library(anndataR)

if(!exists("sample_name")) {
  sample_name <- outFile
}

obj_file=fread(parSampleFile1, header=FALSE)$V1[1]

cat("Reading RDS file:", obj_file, " ...\n")
obj=read_object(obj_file=obj_file, sample_name=sample_name)

cat("Convert to SingleCellExperiment ...\n")
sce <- as.SingleCellExperiment(obj)

cat("Save to h5ad file ...\n")
hd5ad_file <- paste0(sample_name, ".h5ad")
if(file.exists(hd5ad_file)) {
  unlink(hd5ad_file)
  cat("Warning: h5ad file already exists. Overwriting:", hd5ad_file, "\n")
}
write_h5ad(sce, paste0(sample_name, ".h5ad"))

cat("Done.\n")
