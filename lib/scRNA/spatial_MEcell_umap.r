rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/visiumhd/20260206_paula_11498/example_samples/20260209_extract_metadata/spatial_MEcell_umap/result/WHY_01')

### Parameter setting end ###

# Load libraries sf first, otherwise it will cause error when subset seurat object with spatial polygon assay
library(sf)
library(Seurat)
library(logger)
library(arrow)
library(uwot)

log_appender(appender_tee(paste0(sample_name, ".log")))

source("scRNA_func.r")
load_install("MEcell", "liuqivandy/MEcell")

fmap=fread(parSampleFile1, header=FALSE)
obj_file=fmap$V1[1]

cmap=fread(parSampleFile2, header=FALSE)
celltype_column=cmap$V1[1]

k_target=5

log_info(paste0("Read object from ", obj_file, "..."))
obj=readRDS(obj_file)

log_info("Update Seurat object ...")
obj=UpdateSeuratObject(obj)

umap_rds = paste0(sample_name, '.MEcell_umap.k', k_target, '.rds')

mecell=obj$MEcell
while(1){
  row_counts <- Matrix::rowSums(mecell > 0)
  keep_idx <- which(row_counts >= k_target)
  if(length(keep_idx) == nrow(mecell)){
    break
  }
  mecell <- mecell[keep_idx, keep_idx]
}

log_info(paste0("Calculate UMAP with k=", k_target, " ..."))
umap_results <- uwot::umap(
  X = mecell,
  n_neighbors = k_target,
  n_components = 2,
  verbose = TRUE
)

colnames(umap_results) = c("UMAP_1", "UMAP_2")
saveRDS(umap_results, umap_rds)

log_info( paste0(ncol(obj) - nrow(umap_results), " cells with less than ", k_target, " neighbors were removed."))

umap_obj=subset(obj, cells=rownames(umap_results))
stopifnot(all(rownames(umap_results)==colnames(umap_obj)))

umap_obj[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_results),
  key = "UMAP_",
  assay = DefaultAssay(umap_obj)
)

g=get_dim_plot_labelby(umap_obj, reduction='umap', label.by=celltype_column) +
  ggtitle(sample_name)
umap_png=paste0(sample_name, '.MEcell_umap.by_', celltype_column, '.png')
ggsave(umap_png, g, width=10, height=7, dpi=300, units='in', bg='white')

