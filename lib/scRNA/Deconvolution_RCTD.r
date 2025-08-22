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
library(Matrix)
library(tidyverse)
library(hdf5r) # required to read in data file
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)
library(Seurat)
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
spatial_so <- Seurat::Load10X_Spatial(bin.size = 8, data.dir = data_dir, slice = 'slice1')

coords <- Seurat::GetTissueCoordinates(spatial_so)
coords <- coords[,c("x","y")] # keep x and y only

coords$y <- coords$y*-1   # Flip y so that image is not upside down on plots

nUMI <- spatial_so$nCount_Spatial.008um

spatial_counts <- GetAssayData(spatial_so, assay="Spatial.008um", layer="counts")
spatialRNA <- spacexr::SpatialRNA(coords, spatial_counts, nUMI) # 3489

RCTD_obj <- spacexr::create.RCTD(spatialRNA, reference, max_cores = 8)
RCTD_obj <- run.RCTD(RCTD_obj, doublet_mode = "doublet")

saveRDS(RCTD_obj, file = paste0(sample_name, ".RCTD.RDS"))

# process the RCTD results
results <- RCTD_obj@results
# get the cell type proportions
# a data frame of cell type weights for each pixel
cell_type_composition <- results$weights
norm_weights = normalize_weights(results$weights)  # sum as 1 for each spot
cell_type_names <- RCTD_obj@cell_type_info$info[[2]]

# plot cell type composition
# output = cell_type_weights.pdf, each page is a cell type
log_msg("Plotting cell type weights", log_file = log_file)
plot_weights(cell_type_names, spatialRNA, pw_figure, norm_weights)

# Plots all weights for each cell type as in full_mode. (saved as
# output = cell_type_weights_unthreshold.pdf
# only difference is even a spot with very low weight e,g, almost 0, will still be ploted
# plot_weights_unthreshold(cell_type_names, spatialRNA, pw_figure, norm_weights)

# all cell types in a single page
# will use the first_type column in results$results_df as the cell type of a spot
# there are 30 colors, which exceed the default color map (support up to 22 colors), need to add more colors
# and modify the code of the function

n <- length(cell_type_names)
cols <- rainbow(n)
names(cols) <- cell_type_names
results_df <- results$results_df
coords <- spatialRNA@coords

#log_msg("Plotting all cell types", log_file = log_file)
my_plot_all_cell_types <- function(results_df, coords, cell_type_names, resultsdir) {
  barcodes = rownames(results_df[results_df$spot_class != "reject" &
                                   results_df$first_type %in% cell_type_names, ])
  my_table = coords[barcodes, ]
  my_table$class = results_df[barcodes, ]$first_type
  n_levels = length(levels(my_table$class))
  pres = unique(as.integer(my_table$class))
  pres = pres[order(pres)]
  if (n_levels <=21)
    my_pal = pals::kelly(n_levels + 1)[2:(n_levels + 1)]
  if (n_levels > 21)
    my_pal = pals::polychrome(n_levels)
  if (n_levels > 36)
    stop("Plotting currently supports at most 36 cell types as colors")
  
  names(my_pal) <- levels(my_table$class)
  # plot <- ggplot2::ggplot(my_table, ggplot2::aes(x = x, y = y)) +
  #     ggplot2::geom_point(ggplot2::aes(size = 0.15, shape = 19,
  #         color = class)) + ggplot2::scale_color_manual(values = my_pal[pres]) +
  #     ggplot2::scale_shape_identity() + ggplot2::theme_bw() +
  #     ggplot2::scale_size_identity()
  plot = ggplot(my_table, aes(x=x, y=y, color=class)) +
    geom_point(size=0.15, shape=19) +
    scale_color_manual(values = my_pal) +
    theme_bw() +
    theme(aspect.ratio = 1,
          axis.title = element_blank())
}
plot = my_plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, pw_figure)
#increase the legend dot size
ggsave(
  filename = paste0(sample_name, ".all_cell_types.png"),
  plot     = plot,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300,
  bg = "white"
)

# merge the spots, exclude spolts not in RCTD
rctd_spots = rownames(results$weights)
all_spots <- colnames(spatial_so)
common_spots <- intersect(all_spots, rctd_spots)
spatial_so <- subset(spatial_so, cells = common_spots) # 17943  3488

weights_df <- as.data.frame(norm_weights)
colnames(weights_df) <- paste0(colnames(weights_df), "_RCTD")
rctd_features <- colnames(weights_df)
spatial_so <- AddMetaData(spatial_so, metadata = weights_df)
log_msg('Saving spatial_so with RCTD weights', log_file = log_file)
saveRDS(spatial_so, file = paste0(sample_name, ".post_RCTD.RDS"))
