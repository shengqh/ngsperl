rm(list=ls()) 
sample_name='Control_1'
outFile='Control_1'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/shengq2/20260427_14586_HW/cellbender_nd_raw_qc_sct2_h5ad/result/Control_1')

### Parameter setting end ###

source("scRNA_func.r")

library(Seurat)
library(SingleCellExperiment)
load_install("zellkonverter")

if(!exists("sample_name")) {
  sample_name <- outFile
}

obj_file=fread(parSampleFile1, header=FALSE)$V1[1]

cat("Reading RDS file:", obj_file, " ...\n")
obj=read_object(obj_file=obj_file, sample_name=sample_name)

cat("Convert to SingleCellExperiment ...\n")
sce <- as.SingleCellExperiment(obj)

cat("Save to h5ad file ...\n")
writeH5AD(sce, paste0(sample_name, ".h5ad"))

cat("Done.\n")

