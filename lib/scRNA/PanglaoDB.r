rm(list=ls()) 
sample_name='S01_CD3'
outFile='S01_CD3'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20260902_14817_Visium_Complexes_scRNA/raw_qc_PanglaoDB/result/S01_CD3')

### Parameter setting end ###

source("scRNA_func.r")
source("reportFunctions.R")
library(Seurat)
library(SeuratData)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

# possible on dataset level instead of sample level.
if(!exists("sample_name")){
  sample_name=outFile
}

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

bubblemap_width=to_numeric(myoptions$bubblemap_width, 3000)
bubblemap_height=to_numeric(myoptions$bubblemap_height, 1500)
bubblemap_unit=ifelse(bubblemap_width > 50, "px", "in")

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = myoptions$species,
                             curated_markers_file = myoptions$curated_markers_file,
                             HLA_panglao5_file = myoptions$HLA_panglao5_file,
                             layer="Layer4",
                             remove_subtype_str = "",
                             combined_celltype_file = NULL)

tiers = ctdef$tiers

cell_activity_database<-ctdef$cell_activity_database

cat("Cell type annotation by PanglaoDB in cell level ...\n")
data.norm=GetAssayData(obj, assay = "RNA", layer="data")

max_cta_df<-ORA_celltype_cell(cell_exp_data=data.norm,
                              cellType=cell_activity_database$cellType,
                              weight=cell_activity_database$weight)

stopifnot(all(colnames(obj) == max_cta_df$cell))

obj@meta.data$PanglaoDB = max_cta_df$celltype
obj@meta.data$PanglaoDB_cta_score = max_cta_df$cta_score
ct_name="PanglaoDB"

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

major_obj=get_category_with_min_percentage(obj, ct_name, 0.01)
ct_name_count = paste0(ct_name, "_count")
major_obj@meta.data = add_column_count(major_obj@meta.data, ct_name, ct_name_count)

if("umap" %in% names(major_obj@reductions)){
  g=get_dim_plot_labelby(major_obj, label.by = ct_name_count, reduction="umap", pt.size=0.1) + theme(plot.title=element_blank())
  ggsave(paste0(outFile, ".PanglaoDB.qc_umap.png"), g, width=6, height=4, units="in", dpi=300, bg="white")
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
  ggsave(paste0(outFile, ".PanglaoDB.dot.png"), g, width=bubblemap_width, height=bubblemap_height, units=bubblemap_unit, dpi=300, bg="white")
}

rm(major_obj)

if(dir.exists(".local")){
  unlink(".local", recursive=TRUE)
}

