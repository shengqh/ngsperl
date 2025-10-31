rm(list=ls()) 
sample_name='S01_ClassPTC_BRAF'
outFile='S01_ClassPTC_BRAF'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/12904_RB_VisiumHD/20251014_12904_VisiumHD_cellsegment/MEcell_cluster_Leiden/result/S01_ClassPTC_BRAF')

### Parameter setting end ###

source("scRNA_func.r")
source("reportFunctions.R")
library('sf')
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(logger)

log_appender(appender_tee(paste0(sample_name, ".log")))

options(future.globals.maxSize = 600000 * 1024^2)
options(future.seed = TRUE)

myoptions<-read_file_map("fileList2.txt", do_unlist=FALSE)
email=myoptions$email
affiliation=get_affiliation(myoptions)
assay=myoptions$assay
species=myoptions$species

resolution_str=myoptions$MEcell_resolution
resolution=as.numeric(resolution_str)

file_map<-read_file_map("fileList1.txt", do_unlist=FALSE)
sample_name=names(file_map)[1]
rds_file=file_map[[sample_name]]

log_info(paste0("Load Seurat object from ", rds_file, " ..."))
obj = readRDS(rds_file)

stopifnot(assay %in% names(obj@assays))

DefaultAssay(obj) <- assay
log_info(paste0("There are ", ncol(obj), " spots in the Seurat object."))

algorithm=as.numeric(myoptions$cluster_algorithm)
algorithm_name=myoptions$cluster_algorithm_name
log_info(paste0("Clustering algorithm: ", algorithm, "(", algorithm_name, ")"))

log_info(paste0("Perform MEcell clustering with resolution ", resolution, " using algorithm ", algorithm_name))
obj <- FindClusters(obj, resolution = resolution, graph.name="MEcell", algorithm=algorithm, random.seed = 20251026)

res_rds = paste0(sample_name, ".", assay, ".MEcell.res.", resolution, ".", algorithm_name, ".meta.rds")

log_info(paste0("Save Seurat object metadata with MEcell clustering results to ", res_rds))
saveRDS(obj@meta.data, res_rds)

log_info("Done.")
