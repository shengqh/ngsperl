rm(list=ls())
sample_name='Aorta_3019_HIVNeg'
outFile='Aorta_3019_HIVNeg'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/shengq2/20250420_P12891-P12795_10Flex_hg38_cellbender/cellbender_03_clean/result/Aorta_3019_HIVNeg')

### Parameter setting end ###

source("reportFunctions.R")
source("scRNA_func.r")
#We want to keep the cellranger filtered cells, but with cellbender corrected counts.

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
common_genes=intersect(rownames(cellbender_counts), rownames(cellranger_counts))

if(parSampleFile3 != ""){
  decontX_rds=(fread(parSampleFile3, header=FALSE) %>% filter(V2==sample_name))$V1[1]
  cat("Reading ", decontX_rds, " ...\n")
  decontX_obj<-read_scrna_data(decontX_rds)
  decontX_counts=decontX_obj$counts
  cat("decontX_counts: ", nrow(decontX_counts), "genes and", ncol(decontX_counts), "cells.\n")

  common_cells=intersect(common_cells, colnames(decontX_counts))
  df=data.frame(sample=sample_name, 
                n_cell_cellranger=ncol(cellranger_counts), 
                n_cell_cellbender=ncol(cellbender_counts),  
                n_cell_decontX=ncol(decontX_counts),
                n_cell_common=length(common_cells),
                n_gene_cellranger=nrow(cellranger_counts),
                n_gene_cellbender=nrow(cellbender_counts),
                n_gene_common=length(common_genes))
}else{
  df=data.frame(sample=sample_name, 
                n_cell_cellranger=ncol(cellranger_counts), 
                n_cell_cellbender=ncol(cellbender_counts),  
                n_cell_common=length(common_cells),
                n_gene_cellranger=nrow(cellranger_counts),
                n_gene_cellbender=nrow(cellbender_counts),
                n_gene_common=length(common_genes))
}

final_counts=cellbender_counts[common_genes,common_cells]
cat("final counts: ", nrow(final_counts), "genes and", ncol(final_counts), "cells.\n")

write10xCounts( paste0(sample_name, ".cellbender_filtered.clean.h5"), final_counts)

write.table(df, paste0(sample_name, ".clean_summary.txt"), row.names=FALSE, sep="\t")

if(length(names(cellranger_obj)) > 1){
  other_names = names(cellranger_obj)[!(names(cellranger_obj) %in% c("counts"))]
  for(other_name in other_names){
    other_counts = cellranger_obj[[other_name]]
    if(!is.null(other_counts)){
      cat(paste0("cellranger_", other_name), ":", nrow(other_counts), "features and", ncol(other_counts), "cells.\n")
      new_other_counts = other_counts[,common_cells]
      write10xCounts( paste0(sample_name, ".cellbender_filtered.clean.", other_name, ".h5"), new_other_counts)
      cat("final", paste0("cellranger_", other_name), ":", nrow(final_counts), "genes and", ncol(final_counts), "cells.\n")
    }
  }
}
