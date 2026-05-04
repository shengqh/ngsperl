rm(list=ls()) 
outFile='P13667_mix_genome'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/test/20260305_mix_genome_scRNA/mix_genome_count_ratio/result')

### Parameter setting end ###

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

myoptions=fread('fileList2.txt', sep="\t", header=FALSE)
species=myoptions |> dplyr::filter(V2=="species") |> dplyr::pull(V1)

h5_file_list=fread('fileList1.txt', header=FALSE)

all_species_count_per_cell=NULL
ratio_name=""
for(i in 1:nrow(h5_file_list)) {
  if(!file.exists(h5_file_list$V1[i])) {
    stop(paste("File", h5_file_list$V1[i], "does not exist"))
  }

  fname=h5_file_list$V2[i]
  fpath=h5_file_list$V1[i]

  cat("Processing ", fname, ":", fpath, "...\n")
  fdata=Read10X_h5(fpath)
  species_names=sort(unique(gsub("_.+", "", rownames(fdata))))

  gene_species=gsub("_.+", "", rownames(fdata))

  species_count_per_cell=do.call(cbind, lapply(species_names, function(sp){
    Matrix::colSums(fdata[gene_species == sp, , drop=FALSE])
  }))
  colnames(species_count_per_cell)=species_names
  rownames(species_count_per_cell)=colnames(fdata)
  species_count_per_cell=as.data.frame(species_count_per_cell)

  species_count_per_cell$Total=Matrix::rowSums(species_count_per_cell)
  species_count_per_cell$Ratio=species_count_per_cell[, 1] / species_count_per_cell$Total
  species_count_per_cell$Sample=fname

  all_species_count_per_cell=rbind(all_species_count_per_cell, species_count_per_cell)

  ratio_name=paste0(species_names[1], "/TotalReads")
}

g=ggplot(all_species_count_per_cell, aes(x=Ratio)) +
  geom_histogram(bins=30) +
  facet_wrap(~Sample, scale="free_y") +
  labs(title=paste0("Species Count Ratio: ", ratio_name)) +
  theme_bw() +
  theme(strip.background = element_blank(), title=element_text(size=12, face="bold", vjust=0.5, hjust=0.5))

ggsave(paste0(outFile, ".species_count_ratio_histogram.png"), g, width=8, height=6, dpi=300, units="in", bg="white")
