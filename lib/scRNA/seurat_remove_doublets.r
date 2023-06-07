rm(list=ls()) 
sample_name='ENCDO068KYD'
outFile='ENCDO068KYD'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20230504_encode_liver_scRNA_atac_sct2/decontX_nodoublets/result/ENCDO068KYD')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(digest)
library(patchwork)
library(sparseMatrixStats)
library(data.table)
library(tidyr)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

myoptions<-read_file_map(parSampleFile2, do_unlist = FALSE)
doublet_column<-myoptions$doublet_column

#read raw count dat
file_map<-read_file_map(parSampleFile1)
file_path = file_map[[sample_name]]

#read doublet dat
doublet_map<-read_file_map(parSampleFile3)
doublet_file=doublet_map[[sample_name]]

cat("reading", sample_name, ":", file_path, "\n")
lst = read_scrna_data(file_path, keep_seurat_object)

counts = lst$counts  

smeta = readRDS(doublet_file)
dcells = rownames(smeta)[smeta[,doublet_column] %in% c("Doublet", "doublet")]

if(!all(dcells %in% colnames(counts))){
  stop(paste0("Not all doublets cells in file ", file_path))
}

filtered_counts<-counts[,!(colnames(counts) %in% dcells)]

saveRDS(filtered_counts, paste0(sample_name, ".nodoublets.counts.rds"))
