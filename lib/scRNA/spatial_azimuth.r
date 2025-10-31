rm(list=ls()) 
sample_name='S01_ClassPTC_BRAF'
outFile='S01_ClassPTC_BRAF'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/12904_RB_VisiumHD/20251014_12904_VisiumHD_cellsegment/Azimuth/result/S01_ClassPTC_BRAF')

### Parameter setting end ###

source("scRNA_func.r")
source("reportFunctions.R")
load_install("Seurat", "satijalab/seurat")
load_install("Azimuth", "satijalab/azimuth")

options(future.globals.maxSize= 10779361280)
random.seed=20200107

myoptions<-read_file_map(parSampleFile2, do_unlist=FALSE)
assay=myoptions$assay
is_polygons=assay == "Spatial.Polygons"
bin.size=ifelse(is_polygons, "polygons", 8)
assay_slice=ifelse(is_polygons, "slice1.polygons", "slice1.008um")

min_umi=as.numeric(myoptions$nCount_cutoff)

bubblemap_width=to_numeric(myoptions$bubblemap_width, 3000)
bubblemap_height=to_numeric(myoptions$bubblemap_height, 1500)
bubblemap_unit=ifelse(bubblemap_width > 50, "px", "in")

if(file.exists(parSampleFile3)) {
  Azimuth_ref_tbl = fread(parSampleFile3, header=F, stringsAsFactors = F)
  if(sample_name %in% Azimuth_ref_tbl$V2){
    Azimuth_ref_map = split(Azimuth_ref_tbl$V1, Azimuth_ref_tbl$V2)
  }else{
    Azimuth_ref_map = split(Azimuth_ref_tbl$V2, Azimuth_ref_tbl$V1)
  }
  Azimuth_ref = Azimuth_ref_map[[sample_name]]
}else{
  Azimuth_ref=NULL
}

if(is.null(Azimuth_ref)){
  Azimuth_ref=myoptions$Azimuth_ref
}

if(is.null(Azimuth_ref)){
  stop("Azimuth_ref is not defined")
}

cat("Azimuth_ref: ", Azimuth_ref, "\n")

data_dir <- fread(parSampleFile1, header=FALSE)$V1[1]
cache_key = substr(digest::digest(paste0(data_dir, Azimuth_ref), algo="md5"), 1, 8)
cache_rds = paste0(outFile, ".Azimuth.cache.", cache_key, ".rds")
if(file.exists(cache_rds)) {
  cat("Loading cached Azimuth RDS:", cache_rds, "...\n")
  obj <- readRDS(cache_rds)
} else {
  cat("Loading spatial data from:", data_dir, "...\n")
  if(grepl("\\.rds$", tolower(data_dir))) {
    spatial_so <- readRDS(data_dir)
    DefaultAssay(spatial_so) <- assay
  } else {
    spatial_so <- Seurat::Load10X_Spatial(bin.size = bin.size, data.dir = data_dir, slice = 'slice1')
  }

  cat(paste0("Keep the spots with at least ", min_umi, " UMIs ...\n"))
  counts <- GetAssayData(spatial_so, assay=assay, layer="counts")
  umi_counts <- colSums(counts)
  counts <- counts[, umi_counts >= min_umi]

  options(Seurat.object.assay.version = "v3")
  obj <- CreateSeuratObject(counts = counts)
  obj@meta.data$orig.ident <- sample_name

  cat("There are", nrow(obj), "genes and", ncol(obj), "spots\n")

  obj <- RunAzimuth(query = obj, reference = Azimuth_ref)
  saveRDS(obj, cache_rds)
}

cat("meta.data columns: ", paste0(colnames(obj@meta.data), collapse=", "), "\n")

if("predicted.celltype" %in% colnames(obj@meta.data)) {
  obj@meta.data = obj@meta.data |>
    dplyr::mutate(Azimuth_finest=predicted.celltype)
  ct_name = "Azimuth_finest"
}else {
  anno_columns=grep('predicted.+\\d$', colnames(obj@meta.data), value=TRUE)
  azimuth_cols = c()
  idx=1
  for(idx in 1:length(anno_columns)){
    colname = anno_columns[idx]
    newcolname = paste0("Azimuth_l", idx)
    azimuth_cols = c(azimuth_cols, newcolname)
    obj <-AddMetaData(obj, metadata = obj@meta.data[,colname], col.name = newcolname)
  }
  writeLines(azimuth_cols, paste0(outFile, ".Azimuth_cols.txt"))

  if(any(grepl("predicted.+finest", colnames(obj@meta.data)))) {
    anno_columns=grep("predicted.+finest", colnames(obj@meta.data), value=TRUE)
    finest_column=unique(gsub('.score$','',anno_columns))[1]
  }else{
    if("Azimuth_l2" %in% azimuth_cols){
      finest_column = "Azimuth_l2"
    }else if('Azimuth_l1' %in% azimuth_cols) {
      finest_column = "Azimuth_l1"
    }else {
      stop("Cannot find column Azimuth_l1")
    }
  }
  ct_name = "Azimuth_finest"
  obj <-AddMetaData(obj, metadata=obj@meta.data[,finest_column], col.name=ct_name)
}

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

df<-obj@meta.data[,c("orig.ident", ct_name)] %>% 
  dplyr::rename(Sample=orig.ident)
df_tbl<-table(df[,ct_name],df$Sample)
write.csv(df_tbl, paste0(outFile, ".Azimuth_Sample.csv"))

major_obj<-subset(obj, cells=colnames(obj)[!is.na(obj@meta.data[,ct_name])])
rm(obj)

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

major_obj=get_category_with_min_percentage(major_obj, ct_name, 0.01)
ct_name_count = paste0(ct_name, "_count")
major_obj@meta.data = add_column_count(major_obj@meta.data, ct_name, ct_name_count)

g=get_dim_plot_labelby(major_obj, label.by = ct_name, reduction="ref.umap", pt.size=0.1) + theme(plot.title=element_blank())
ggsave(paste0(outFile, ".Azimuth.png"), g, width=6, height=4, units="in", dpi=300, bg="white")

if("umap" %in% names(major_obj@reductions)){
  g=get_dim_plot_labelby(major_obj, label.by = ct_name_count, reduction="umap", pt.size=0.1) + theme(plot.title=element_blank())
  ggsave(paste0(outFile, ".Azimuth.qc_umap.png"), g, width=6, height=4, units="in", dpi=300, bg="white")
}

if(has_bubblemap){
  g<-get_bubble_plot(
    obj=major_obj, 
    cur_res=NA, 
    cur_celltype=ct_name_count, 
    bubblemap_file, 
    assay="RNA", 
    species=myoptions$species,
    dot.scale=4)
  ggsave(paste0(outFile, ".Azimuth.dot.png"), g, width=bubblemap_width, height=bubblemap_height, units=bubblemap_unit, dpi=300, bg="white")
}

rm(major_obj)

if(dir.exists(".local")){
  unlink(".local", recursive=TRUE)
}
