rm(list=ls()) 
sample_name='CW_1'
outFile='CW_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/celly_wanjalla_projects/20260209_14382_scRNA_cellranger/20260226_mix_genome_scRNA/split_mix_genome_counts_hg38/result/CW_1')

### Parameter setting end ###

library(Seurat)
library(data.table)
library(dplyr)

myoptions=fread('fileList2.txt', header=FALSE)
species=myoptions |> dplyr::filter(V2=="species") |> dplyr::pull(V1)

barcodes_file=fread('fileList1.txt', header=FALSE)$V1[1]
barcodes=fread(barcodes_file, header=FALSE) |>
  dplyr::rename(Species = V1, Cell=V2)

if(!(species %in% barcodes$Species)) {
  stop(paste("Species", species, "not found in barcodes"))
}
cells <- barcodes |> dplyr::filter(Species == species) |> dplyr::pull(Cell)

species_h5_file = fread('fileList3.txt', header=FALSE) |>
  dplyr::filter(V2 == sample_name) |> dplyr::pull(V1)

counts=Read10X_h5(species_h5_file)

filtered_counts <- counts[, colnames(counts) %in% cells]

cat(ncol(filtered_counts), "out of", length(cells), "cells successfully extracted for species", species, "\n")

saveRDS(filtered_counts, file = paste0(sample_name, ".", species, ".counts.rds"))
