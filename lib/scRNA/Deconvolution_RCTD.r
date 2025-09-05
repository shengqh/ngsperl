rm(list=ls()) 
sample_name='S03_ClassPTC_BRAF'
outFile='S03_ClassPTC_BRAF'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/h_cqs/jstolze/VWeiss/20250806_Visium_RCTD_Deconvolute/result/20250807_FastMNN_RCTD_Reference_general.RDS'
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20250812_VisiumHD_RCTD_test/RCTD/result/S03_ClassPTC_BRAF')

### Parameter setting end ###

source("Deconvolution_functions.R")
library(Seurat)
library(dplyr)
library(spacexr)
library(hdf5r) # required to read in data file
library(data.table)

options_df=fread("fileList2.txt", header=FALSE)
myoptions=split(options_df$V1, options_df$V2)

#use 8um right now.
#myoptions$bin.size=as.numeric(myoptions$bin.size)
myoptions$RCTD_thread=as.numeric(myoptions$RCTD_thread)

log_file = paste0(sample_name, ".RCTD.log")
ignored = unlink(log_file, force=TRUE)

log_msg(paste0("Starting RCTD deconvolution for sample: ", sample_name), log_file = log_file)

log_msg(paste0("Loading reference RDS file: ", parFile1), log_file = log_file)
reference=readRDS(parFile1)

data_dir <- fread(parSampleFile1, header=FALSE)$V1[1]
log_msg(paste0("Loading spatial data from: ", data_dir), log_file = log_file)

if(grepl("\\.rds$", tolower(data_dir))) {
  spatial_so <- readRDS(data_dir)
} else {
  spatial_so <- Seurat::Load10X_Spatial(bin.size = 8, data.dir = data_dir, slice = 'slice1')
}

coords <- Seurat::GetTissueCoordinates(spatial_so)
coords <- coords[,c("x","y")] # keep x and y only

nUMI <- spatial_so$nCount_Spatial.008um

spatial_counts <- GetAssayData(spatial_so, assay="Spatial.008um", layer="counts")
spatialRNA <- spacexr::SpatialRNA(coords, spatial_counts, nUMI)

rctd_rds = paste0(sample_name, ".RCTD.RDS")
if(!file.exists(rctd_rds)){
  RCTD_obj <- spacexr::create.RCTD(spatialRNA, reference, max_cores = myoptions$RCTD_thread)
  RCTD_obj <- run.RCTD(RCTD_obj, doublet_mode = "doublet")
  saveRDS(RCTD_obj, file = rctd_rds)
}else{
  RCTD_obj <- readRDS(rctd_rds)
}

# process the RCTD results
results <- RCTD_obj@results

# get the cell type proportions
# a data frame of cell type weights for each pixel
# sum to 1, recommended from RCTD author, https://github.com/dmcable/spacexr/issues/45
weights_df = normalize_weights(results$weights) |> as.data.frame()  # sum as 1 for each spot
colnames(weights_df) <- paste0("RCTD2_", colnames(weights_df))

# merge the spots, exclude spolts not in RCTD
common_spots <- intersect(colnames(spatial_so), rownames(weights_df))
spatial_so <- subset(spatial_so, cells = common_spots)

# Add weights into spatial_so meta data
spatial_so <- AddMetaData(spatial_so, metadata = weights_df)

# Add results from doublet mode into spatial_so meta data
rctd_df=results$results_df |> dplyr::select("spot_class", "first_type", "second_type")
names(rctd_df) <- c("RCTD1_Class", "RCTD1_Label1", "RCTD1_Label2")
spatial_so <- AddMetaData(spatial_so, metadata = rctd_df)

log_msg('Saving spatial_so with RCTD weights', log_file = log_file)
saveRDS(spatial_so, file = paste0(sample_name, ".post_RCTD.RDS"))
