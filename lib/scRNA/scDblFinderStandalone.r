library(SingleCellExperiment)
library(Seurat)
library(scDblFinder)

args <- commandArgs(TRUE)

if(length(args) == 0){
  h5file = "/data/vickers_lab/2022/20221201_9061_DM/Count/9061-DM-1/filtered_feature_bc_matrix.h5"
  outFile = "DM_1"
}else{
  h5file = args[1]
  outFile = args[2]
}

counts <- Read10X_h5(h5file)
if(is.list(counts)){
  counts<-counts$`Gene Expression`
}

sce <- SingleCellExperiment(list(counts = counts))
rm(counts)

set.seed(20230118)
sce <- scDblFinder(sce)

obj <- as.Seurat(sce, data = NULL)
rm(sce)

obj$orig.ident<-outFile

saveRDS(obj@meta.data, paste0(outFile, ".scDblFinder_meta.rds"))

write.csv(table(obj$scDblFinder.class), paste0(outFile, ".scDblFinder.csv"))

singlet <- subset(obj, subset = scDblFinder.class == "singlet")
saveRDS(singlet, paste0(outFile, ".scDblFinder.rds"))

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')

sce<-as.SingleCellExperiment(obj)
plotDoubletMap(sce)
