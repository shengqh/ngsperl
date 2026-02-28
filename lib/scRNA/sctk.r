rm(list=ls()) 
sample_name='Adipose_9240'
outFile='Adipose_9240'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230410_combined_scRNA_hg38_test_vst2/decontX_sctk/result/Adipose_9240')

### Parameter setting end ###

source("scRNA_func.r")
library("singleCellTK")
library("SingleCellExperiment")
library("Seurat")

sample_map = read_file_map("fileList1.txt")

myoptions = read_file_map("fileList2.txt", do_unlist=FALSE)
nFeature_cutoff_min=max(as.numeric(myoptions$nFeature_cutoff_min) - 50, 100)
nCount_cutoff=max(as.numeric(myoptions$nCount_cutoff) - 50, 200)

cat("sample_name\n")
countfile = sample_map[sample_name]

cat("  read_scrna_data\n")
lst = read_scrna_data(countfile)
counts<-lst$counts
rawobj = CreateSeuratObject(counts = counts, project = sample_name)
rawobj<-subset(rawobj, subset = nFeature_RNA >= nFeature_cutoff_min & nCount_RNA >= nCount_cutoff)

sce <- as.SingleCellExperiment(rawobj)
rm(rawobj)

cat("  runCellQC\n")
#sce = runCellQC(sce, algorithms= c("QCMetrics", "scrublet", "doubletFinder", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid", "decontX"))
#There is error to run scrublet, so ignore it.
sce = runCellQC(sce, algorithms= c("QCMetrics", "doubletFinder", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid", "decontX"))

saveRDS(colData(sce), paste0(sample_name, ".meta.rds"))
