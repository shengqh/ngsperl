rm(list=ls()) 
sample_name='CD_Met_01'
outFile='CD_Met_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/maureen_gannon_projects/20240320_10940_snRNAseq_mmulatta_proteincoding_cellbender/cellbender_03_clean/result/CD_Met_01')

### Parameter setting end ###

#We want to keep the cellranger filtered cells, but with cellbender corrected counts.

source("scRNA_func.r")
library(DropletUtils)

cellbender_filtered_h5=fread(parSampleFile1, header=FALSE)$V1[1]
cat("Reading ", cellbender_filtered_h5, " ...\n")
cellbender_obj<-read_scrna_data(cellbender_filtered_h5)
cellbender_counts=cellbender_obj$counts
cat("cellbender_counts: ", nrow(cellbender_counts), "genes and", ncol(cellbender_counts), "cells.\n")

cellranger_filtered_h5=(fread(parSampleFile2, header=FALSE) %>% filter(V2==sample_name))$V1[1]
cat("Reading ", cellranger_filtered_h5, " ...\n")
cellranger_obj<-read_scrna_data(cellranger_filtered_h5)
cellranger_counts=cellranger_obj$counts
cat("cellranger_counts: ", nrow(cellranger_counts), "genes and", ncol(cellranger_counts), "cells.\n")

common_cells=intersect(colnames(cellbender_counts), colnames(cellranger_counts))

final_counts=cellbender_counts[,common_cells]
cat("final counts: ", nrow(final_counts), "genes and", ncol(final_counts), "cells.\n")

write10xCounts( paste0(sample_name, ".cellbender_filtered.clean.h5"), final_counts)

df=data.frame(sample=sample_name, 
              n_cellranger=ncol(cellranger_counts), 
              n_cellbender=ncol(cellbender_counts),  
              n_common=length(common_cells))
write.table(df, paste0(sample_name, ".clean_summary.txt"), row.names=FALSE, sep="\t")