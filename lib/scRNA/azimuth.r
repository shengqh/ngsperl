rm(list=ls()) 
sample_name='cvd_11a'
outFile='cvd_11a'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250211_T04_snRNA_hg38/raw_qc_Azimuth/result/cvd_11a')

### Parameter setting end ###

source("reportFunctions.R")
source("scRNA_func.r")
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

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

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

if(myoptions$species == "Mm") {
  if(require('homologene')) {
    cat("Source data is mouse, convert to human gene symbols\n")
    mm2hs = mouse_symbol_to_human_symbol(rownames(obj))

    # Rename and (optionally) collapse duplicates
    common <- intersect(rownames(obj), names(mm2hs))

    cat("Common genes: ", length(common), "\n")

    # --- 0) Prep
    library(Seurat)
    library(Matrix)

    mat <- GetAssayData(obj, layer = "counts")

    # Determine which genes can be renamed
    common <- intersect(rownames(mat), names(mm2hs))
    mat <- mat[common, ]
    rownames(mat) = mm2hs[common]

    new_obj <- CreateSeuratObject(counts = mat)
    # must add one more column to meta.data, otherwise, scTransform will fail.
    new_obj <- AddMetaData(new_obj, metadata = "azimuth", col.name = "project")
    new_obj <- SCTransform(new_obj, verbose = FALSE)

    new_obj <- RunAzimuth(query = new_obj, reference = Azimuth_ref)
    obj@meta.data <- cbind(obj@meta.data, new_obj@meta.data |> dplyr::select(-orig.ident, -project))
    obj[["ref.umap"]] <- new_obj[["ref.umap"]]
  }else{
    cat("homologene is not installed, will not convert mouse symbols to human gene symbols, can use RNA assay only.\n")
    obj <- RunAzimuth(query = obj, reference = Azimuth_ref, assay="RNA")
  }
} else {
  obj <- RunAzimuth(query = obj, reference = Azimuth_ref)
}

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
  }else{
    finest_column = "Azimuth_l1"
  }
}
ct_name = "Azimuth_finest"
obj <-AddMetaData(obj, metadata=obj@meta.data[,finest_column], col.name=ct_name)

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

