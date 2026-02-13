rm(list=ls()) 
sample_name='WHY_04'
outFile='WHY_04'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/visiumhd/20260206_paula_11498/example_samples/spatial_MEcell_umap_k10/result/WHY_04')

### Parameter setting end ###

source("scRNA_func.r")
# Load libraries sf first, otherwise it will cause error when subset seurat object with spatial polygon assay
library(sf)
library(Seurat)
library(logger)
library(arrow)
library(uwot)

log_appender(appender_tee(paste0(sample_name, ".log")))

load_install("MEcell", "liuqivandy/MEcell")

fmap=fread(parSampleFile1, header=FALSE)
obj_file=fmap$V1[1]

cmap=fread(parSampleFile2, header=FALSE)
celltype_column=cmap |> 
  dplyr::filter(V2==sample_name) |>
  dplyr::pull(V1)

omap=fread(parSampleFile3, header=FALSE)
k_neighbors=omap |> 
  dplyr::filter(V2=="k_neighbors") |>
  dplyr::pull(V1) |>
  as.numeric()

log_info(paste0("Read object from ", obj_file, "..."))
obj=readRDS(obj_file)

if(!celltype_column %in% colnames(obj@meta.data)){
  stop(paste0("Cell type column ", celltype_column, " not found in object meta data!"))
}

log_info("Update Seurat object ...")
obj=UpdateSeuratObject(obj)

mecell=obj$MEcell

# Calculate UMAP using MEcell as graph with all cells
umap_rds = paste0(sample_name, '.MEcell_umap.graph_all.rds')
if(!file.exists(umap_rds)){
  log_info(paste0("Calculate UMAP using MEcell as graph with all cells ..."))
  umap_results <- uwot::umap(
    X = NULL,
    nn_method = mecell,
    n_threads = 8,
    n_sgd_threads = 1, # result won't be perfectly reproducible if n_sgd_threads > 1
    verbose = TRUE
  )

  colnames(umap_results) = c("UMAP_1", "UMAP_2")
  saveRDS(umap_results, umap_rds)
}else {
  umap_results = readRDS(umap_rds)
}

stopifnot(all(rownames(umap_results)==colnames(obj)))

obj[["umapgraphall"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_results),
  key = "UMAP_",
  assay = DefaultAssay(obj)
)

g=get_dim_plot_labelby(obj, reduction='umapgraphall', label.by=celltype_column) +
  ggtitle(sample_name)
umap_png=paste0(sample_name, '.MEcell_umap.graph_all.png')
ggsave(umap_png, g, width=10, height=7, dpi=300, units='in', bg='white')

# Filter cells with less than k_neighbors neighbors
log_info(paste0("Filter cells with less than ", k_neighbors, " neighbors ..."))

mecell=obj$MEcell
while(1){
  row_counts <- Matrix::rowSums(mecell > 0)
  keep_idx <- which(row_counts >= k_neighbors)
  if(length(keep_idx) == nrow(mecell)){
    break
  }
  mecell <- mecell[keep_idx, keep_idx]
}

log_info( paste0(ncol(obj) - nrow(mecell), " cells with less than ", k_neighbors, " neighbors were removed."))
umap_obj=subset(obj, cells=rownames(mecell))

# Calculate UMAP using MEcell as graph after k_neighbors neighbors filtering
umap_rds = paste0(sample_name, '.MEcell_umap.graph_k', k_neighbors, '.rds')
if(!file.exists(umap_rds)){
  log_info(paste0("Calculate UMAP with MEcell as graph ..."))
  umap_results <- uwot::umap(
    X = NULL,
    nn_method = mecell,
    n_threads = 8,
    n_sgd_threads = 1, # result won't be perfectly reproducible if n_sgd_threads > 1
    verbose = TRUE
  )

  colnames(umap_results) = c("UMAP_1", "UMAP_2")

  saveRDS(umap_results, umap_rds)
}else {
  umap_results = readRDS(umap_rds)
}
stopifnot(all(rownames(umap_results)==colnames(umap_obj)))

umap_obj[["umapgraphk"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_results),
  key = "UMAP_",
  assay = DefaultAssay(umap_obj)
)

g=get_dim_plot_labelby(umap_obj, reduction='umapgraphk', label.by=celltype_column) +
  ggtitle(sample_name)
umap_png=paste0(sample_name, '.MEcell_umap.graph_k', k_neighbors, '.png')
ggsave(umap_png, g, width=10, height=7, dpi=300, units='in', bg='white')

# Calculate UMAP using MEcell as distance with k_neighbors neighbors
umap_rds = paste0(sample_name, '.MEcell_umap.dist_k', k_neighbors, '.rds')
if(!file.exists(umap_rds)){
  log_info(paste0("Calculate UMAP using MEcell as distance with ", k_neighbors, " neighbors ..."))
  umap_results <- uwot::umap(
    X = mecell,
    n_neighbors = k_neighbors,
    n_components = 2,
    n_threads = 8,
    n_sgd_threads = 1, # result won't be perfectly reproducible if n_sgd_threads > 1
    verbose = TRUE
  )

  colnames(umap_results) = c("UMAP_1", "UMAP_2")
  saveRDS(umap_results, umap_rds)
}else {
  umap_results = readRDS(umap_rds)
}

stopifnot(all(rownames(umap_results)==colnames(umap_obj)))
umap_obj[["umapdist"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_results),
  key = "UMAP_",
  assay = DefaultAssay(umap_obj)
)

g=get_dim_plot_labelby(umap_obj, reduction='umapdist', label.by=celltype_column) +
  ggtitle(sample_name)
umap_png=paste0(sample_name, '.MEcell_umap.dist_k', k_neighbors, '.png')
ggsave(umap_png, g, width=10, height=7, dpi=300, units='in', bg='white')


