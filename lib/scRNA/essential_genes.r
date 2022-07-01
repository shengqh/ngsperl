rm(list=ls()) 
outFile='AG3669'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('C:/projects/data/h_gelbard_lab/projects/20220609_scRNA_3669_sct/essential_genes/result')


### Parameter setting end ###

source("scRNA_func.r")


library(Seurat)
library(ggplot2)
library(data.table)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

species=myoptions$species
markerfile<-myoptions$db_markers_file
HLA_panglao5_file=myoptions$HLA_panglao5_file
remove_subtype_of=myoptions$remove_subtype
cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file, curated_markers_file=myoptions$curated_markers_file)
genes<-unlist(cell_activity_database$cellType)

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)
if(has_bubblemap){
  genes_df <- read_bubble_genes(bubblemap_file, NA)
  bubble_genes<-unique(genes_df$gene)
  genes<-c(genes, bubble_genes)
}

if(file.size(parSampleFile2) > 0L){
  marker_gene_files = read.table(parSampleFile2, sep="\t", header=F)$V1
  for (mfile in marker_gene_files){
    mg=read.table(mfile, sep="\t")[,2]
    genes<-c(genes, mg)
  }
}

genes<-unique(genes)
genes<-genes[order(genes)]

writeLines(genes, paste0(outFile, ".txt"))
