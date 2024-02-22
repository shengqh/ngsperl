rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/00_prepare_gex_h5/result')

### Parameter setting end ###

source("scRNA_func.r")
library(DropletUtils)
#for cellbender, if the protein capture data is included, it will be used as genes. so we need to keep gene expression data only

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

h5_map = read_file_map(parSampleFile1)

sample = names(h5_map)[1]
for(sample in names(h5_map)){
  sfile = h5_map[[sample]]
  cat("  read", sfile, "\n")
  obj = Read10X_h5(sfile)
  gex = obj$`Gene Expression`

  write10xCounts(paste0(sample, myoptions$suffix), gex)
}

