rm(list=ls()) 
outFile='P9061'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_rawdata/result/P9061.rawobj.rds'
parFile2=''
parFile3=''


setwd('/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_rawdata_scDblFinder/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(scDblFinder)
library(plyr)

obj <- read_object(parFile1)
sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce, samples = "sample")

obj <- as.Seurat(sce)
rm(sce)

saveRDS(obj@meta.data, paste0(outFile, ".scDblFinder_meta.rds"))

write.csv(table(obj$scDblFinder.class, obj$scDblFinder.sample), paste0(outFile, ".scDblFinder.csv"))

singlet <- subset(obj, subset = scDblFinder.class == "singlet")
saveRDS(singlet, paste0(outFile, ".scDblFinder.rds"))
