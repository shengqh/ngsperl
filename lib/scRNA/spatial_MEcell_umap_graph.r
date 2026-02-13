rm(list=ls()) 
sample_name='S01_ClassPTC_BRAF'
outFile='S01_ClassPTC_BRAF'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/12904_RB_VisiumHD/20260212_12904_VisiumHD_cellsegment/MEcell_umap/result/S01_ClassPTC_BRAF')

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
if("celltype_column" %in% cmap$V2){
  celltype_column=cmap |> dplyr::filter(V2 == "celltype_column") |> dplyr::pull(V1)
} else {
  celltype_column=NULL
}
min_neighbors=as.numeric(cmap |> dplyr::filter(V2 == "min_neighbors") |> dplyr::pull(V1))
min_neighbors=1

log_info(paste0("Read object from ", obj_file, "..."))
obj=readRDS(obj_file)

log_info("Update Seurat object ...")
obj=UpdateSeuratObject(obj)

mecell=obj$MEcell
n_inf=sum(is.infinite(mecell@x))
if(n_inf > 0){
  log_warn(paste0("MEcell matrix contains ", n_inf, " infinite values, which will be set to 0 for UMAP calculation."))
  mecell@x[!is.finite(mecell@x)] <- 0
}
while(1){
  row_counts <- Matrix::rowSums(mecell > 0)
  keep_idx <- which(row_counts >= min_neighbors)
  if(length(keep_idx) == nrow(mecell)){
    break
  }
  mecell <- mecell[keep_idx, keep_idx]
}
n_deleted=ncol(obj) - nrow(mecell)
log_info( paste0(n_deleted, " cells with less than ", min_neighbors, " neighbors were removed."))

log_info(paste0("Calculate UMAP ..."))
set.seed(20260210)
umap_results <- uwot::umap(
  X = NULL,
  n_threads = 8,
  n_sgd_threads = 1, # result won't be perfectly reproducible if n_sgd_threads > 1
  nn_method = mecell,
  verbose = TRUE
)

colnames(umap_results) = c("UMAP_1", "UMAP_2")

umap_rds = paste0(sample_name, '.MEcell_umap.rds')
saveRDS(umap_results, umap_rds)

if(!is.null(celltype_column)){
  if(n_deleted > 0){
    obj = subset(obj, cells=rownames(umap_results))
    stopifnot(all(colnames(obj) == rownames(umap_results)))
  }
  obj[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap_results),
    key = "UMAP_",
    assay = DefaultAssay(obj)
  )

  g=get_dim_plot_labelby(obj, reduction='umap', label.by=celltype_column, legend.title="Cell type") +
    ggtitle(sample_name)

  umap_png=paste0(sample_name, '.MEcell_umap.png')
  ggsave(umap_png, g, width=10, height=7, dpi=300, units='in', bg='white')
}

