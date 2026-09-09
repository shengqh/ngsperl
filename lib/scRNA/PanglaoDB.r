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
                             layer = "Layer4",
                             remove_subtype_str = "",
                             combined_celltype_file = NULL)

tiers = ctdef$tiers

cell_activity_database<-ctdef$cell_activity_database

cat("Cell type annotation by PanglaoDB in cell level ...\n")
data.norm=GetAssayData(obj, assay = "RNA", layer="data")

max_cta_df<-CTA_celltype_cell(cell_exp_data=data.norm,
                              cellType=cell_activity_database$cellType,
                              weight=cell_activity_database$weight)

stopifnot(all(colnames(obj) == max_cta_df$cell))
meta = obj@meta.data

meta$PanglaoDB_cta_raw = max_cta_df$celltype
meta$PanglaoDB_cta_score = max_cta_df$cta_score

meta=meta |> 
  tibble::rownames_to_column("Cell_barcode") |>
  dplyr::left_join(tiers, by=c("PanglaoDB_cta_raw"="Celltype.name")) |>
  tibble::column_to_rownames("Cell_barcode")

layer1=c("Epithelial cells", "Neural cell", "Muscle cell")
meta = meta |>
  dplyr::mutate(PanglaoDB = ifelse(Layer1 %in% layer1, Layer1, Layer2))

ct_name="PanglaoDB"

obj@meta.data=meta

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

marker_genes = cell_activity_database$cellType[unique(major_obj@meta.data$PanglaoDB_cta_raw)]
marker_gene_df = df <- stack(marker_genes) |>
  dplyr::rename("gene"="values", "celltype"="ind") |>
  dplyr::left_join(tiers, by=c("celltype"="Celltype.name")) |>
  dplyr::select(Layer2, gene) |>
  dplyr::distinct() |>
  dplyr::filter(gene %in% row.names(obj))

gene_groups=split(marker_gene_df$gene, marker_gene_df$Layer2)
g=get_dot_plot(obj=major_obj, group.by=ct_name, gene_groups=gene_groups)

# keep the gene with pct.exp > 5% and avg.exp.scaled > 0.1
# then for each gene, keep the top cell type
# then for each cell type, keep the top 5 genes
gdata = g$data |>
  dplyr::mutate(id=as.character(id), feature.groups=as.character(feature.groups)) |>
  dplyr::filter(id == feature.groups) |>
  dplyr::filter(pct.exp > 5) |>
  dplyr::filter(avg.exp.scaled > 0.1) |>
  dplyr::select(id, features.plot, avg.exp) |>
  dplyr::group_by(features.plot) |>  
  dplyr::top_n(1, avg.exp) |>
  dplyr::ungroup() |>
  dplyr::group_by(id) |>
  dplyr::top_n(5, avg.exp) |>
  dplyr::ungroup()

gene_groups=split(gdata$features.plot, gdata$id)
saveRDS(gene_groups, paste0(outFile, ".PanglaoDB.markers.rds"))

g=get_dot_plot( major_obj, 
                "PanglaoDB", 
                gene_groups, 
                assay="RNA", 
                rotate.title=TRUE, 
                use_blue_yellow_red=TRUE, 
                dot.scale=4)

dot_height=get_dot_height_num(length(gene_groups))
dot_width=get_dot_width(g)
ggsave(paste0(outFile, ".PanglaoDB.markers.dot.png"), g, width=dot_width, height=dot_height, units="px", dpi=300, bg="white")

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

