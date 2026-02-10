rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/visiumhd/20260206_paula_11498/example_samples/20260209_extract_metadata/extract_metadata/result/WHY_01')

### Parameter setting end ###

library(Seurat)

fmap=read.table(parSampleFile1, header=F, stringsAsFactors=F, sep="\t")
obj_file=fmap$V1[1]

obj=readRDS(obj_file)
meta=obj@meta.data

saveRDS(meta, paste0(outFile, '.meta.rds'))
writeLines(colnames(meta), paste0(outFile, '.meta_colnames.txt'))
