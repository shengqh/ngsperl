rm(list=ls()) 
sample_name='S4'
outFile='S4'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/pipelines/20260213_cosmx/result/MEcell/result/S4')

### Parameter setting end ###

source("scRNA_func.r")
# Load libraries sf first, otherwise it will cause error when subset seurat object with spatial polygon assay
library(sf)
library(logger)

log_appender(appender_tee(paste0(sample_name, ".log")))

load_install("MEcell", "liuqivandy/MEcell")

myoptions_tbl=fread(parSampleFile2, header=FALSE) 
myoptions=split(myoptions_tbl$V1, myoptions_tbl$V2)

assay=myoptions$assay

data_dir <- fread(parSampleFile1, header=FALSE)$V1[1]

log_info(paste0("Loading spatial data from:", data_dir, "...\n"))
if(grepl("\\.rds$", tolower(data_dir))) {
  log_info(paste0("Reading RDS file:", data_dir, "...\n"))
  object <- readRDS(data_dir)
  DefaultAssay(object) <- assay
} else {
  log_info(paste0("Loading 10X Spatial data from: ", data_dir, "...\n"))
  object <- Seurat::Load10X_Spatial(bin.size = bin.size, data.dir = data_dir, slice = 'slice1')
}

log_info(paste0("There are ", nrow(object), " genes and ", ncol(object), " spots\n"))

log_info("Normalizing data")
object <- NormalizeData(object)

log_info("FindVariableFeatures")
object <- FindVariableFeatures(object)

log_info("ScaleData")
object <- ScaleData(object)

log_info("RunPCA")
object <- RunPCA(object)

log_info("Run MEcell")
object <- MEcell(object, usepca=TRUE)

log_info(paste0("Save Seurat object with MEcell results to ", outFile, ".MEcell.rds"))
saveRDS(object, file=paste0(outFile, ".MEcell.rds"))

mtx_file = paste0(outFile, ".MEcell.mtx")
log_info(paste0("Save MEcell to ", mtx_file, ".gz ..."))
mecell = object$MEcell
sparse_mtx = as(mecell, "dgCMatrix")
ignored = writeMM(obj = sparse_mtx, file = mtx_file)
system(paste0("gzip -f ", mtx_file))

writeLines(rownames(mecell), paste0(outFile, ".MEcell.cells.txt"))

log_info("Done.")

