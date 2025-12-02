rm(list=ls()) 
sample_name='DM_1'
outFile='DM_1'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20251128_9061_scRNA_mm10_cellbender_default_redo/cellbender_decontX_scDblFinder/result/DM_1')

### Parameter setting end ###

source("scRNA_func.r")
source("countTableVisFunctions.R")
library(Seurat)
library(scDblFinder)

h5_file=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)$V1[1]

obj <- read_object(h5_file)

if(file.exists(parSampleFile2)){
  meta_file = read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F) |>
    dplyr::filter(V2==sample_name) |>
    dplyr::pull(V1)

  meta.data = readRDS(meta_file)
  obj = subset(obj, cells=rownames(meta.data))
  stopifnot(all(colnames(obj) == rownames(meta.data)))

  obj@meta.data = meta.data

  clusters = fread(parSampleFile3, header=F, stringsAsFactors = F) |>
    dplyr::filter(V2=="ct_column") |>
    dplyr::pull(V1)
  cat("clusters=", clusters, "\n")
}else{
  clusters = NULL
}

sce <- as.SingleCellExperiment(obj)

sce <- scDblFinder(sce, clusters=clusters)

obj <- as.Seurat(sce, data = NULL)
rm(sce)

saveRDS(obj@meta.data, paste0(outFile, ".scDblFinder.meta.rds"))

obj$orig.ident = sample_name
saveRDS(obj, paste0(outFile, ".scDblFinder.object.rds"))
write.csv(table(obj$scDblFinder.class, obj$orig.ident), paste0(outFile, ".scDblFinder.class.csv"))

singlet <- subset(obj, subset = scDblFinder.class == "singlet")
saveRDS(singlet, paste0(outFile, ".scDblFinder.singlet_object.rds"))

