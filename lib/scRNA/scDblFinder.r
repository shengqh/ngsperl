rm(list=ls())
sample_name='KA_0001'
outFile='KA_0001'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20241030_T01_cellbender/scDblFinder/result/KA_0001')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(scDblFinder)
library(plyr)

h5_file=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)$V1[1]

obj <- read_object(h5_file)
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)

obj <- as.Seurat(sce, data = NULL)
rm(sce)

saveRDS(obj@meta.data, paste0(outFile, ".scDblFinder.meta.rds"))

obj$orig.ident = sample_name
saveRDS(obj, paste0(outFile, ".scDblFinder.object.rds"))
write.csv(table(obj$scDblFinder.class, obj$orig.ident), paste0(outFile, ".scDblFinder.class.csv"))

singlet <- subset(obj, subset = scDblFinder.class == "singlet")
saveRDS(singlet, paste0(outFile, ".scDblFinder.singlet_object.rds"))
