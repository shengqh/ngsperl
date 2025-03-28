rm(list=ls()) 
sample_name='KA_0001'
outFile='KA_0001'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250325_T01_cellbender/cellbender_00_extract_gene_expression_h5/result/KA_0001')

### Parameter setting end ###

source("reportFunctions.R")
library(Seurat)
library(DropletUtils)
#for cellbender, if the protein capture data is included, it will be used as genes. so we need to keep gene expression data only

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

h5_map = read_file_map(parSampleFile1)
sfile = h5_map[[sample_name]]

cat("  read", sfile, "\n")
obj = Read10X_h5(sfile)
gex = obj$`Gene Expression`

write10xCounts(paste0(sample_name, myoptions$suffix), gex)

