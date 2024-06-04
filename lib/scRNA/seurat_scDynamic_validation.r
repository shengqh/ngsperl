rm(list=ls()) 
outFile='coronary'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parFile1='/nobackup/shah_lab/shengq2/20240208_CAC_proteomics_scRNA/chiara_scRNA/20240531_T01_prepare_data_coronary/coronary.DE.rds'
parFile2='/nobackup/shah_lab/shengq2/20240208_CAC_proteomics_scRNA/chiara_scRNA/20240531_T01_prepare_data_coronary/coronary.meta.rds'
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240208_CAC_proteomics_scRNA/chiara_scRNA/20240601_T02_subcluster/celltype_validation/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)

options(future.globals.maxSize= 10779361280)

myoptions = read_file_map(parSampleFile1, do_unlist = FALSE)
doublet_column = myoptions$doublet_column
celltype_column = myoptions$celltype_column
bubblemap_file = myoptions$bubblemap_file
cur_assay = ifelse(is_one(myoptions$by_sctransform), "SCT", "RNA")
species = myoptions$species
create_clusters = is_one(myoptions$create_clusters)

file_dir=paste0(outFile, gsub(".html", "", myoptions$rmd_ext))
dir.create(file_dir, showWarnings=FALSE)
file_prefix=file.path(file_dir, outFile)

meta<-readRDS(parFile2)
if(create_clusters){
  celltype_cluster_column = paste0(celltype_column, "_clusters")
  meta[,celltype_column]=factor_by_count(meta[,celltype_column])
  meta[,celltype_cluster_column]=as.numeric(meta[,celltype_column])-1
}else{
  if(celltype_column == "seurat_cell_type"){
    celltype_cluster_column = "seurat_clusters"
  }else if (paste0("seurat_", celltype_column) %in% colnames(meta)){
    celltype_cluster_column = "seurat_clusters"
  }else{
    celltype_cluster_column = paste0(celltype_column, "_clusters")
  }
}

stopifnot(celltype_cluster_column %in% colnames(meta))

#meta[,celltype_cluster_column] <- as.character(meta[,celltype_cluster_column])

validation_columns=c()

has_decontX = exists('parSampleFile6')
if(has_decontX){
  meta = fill_meta_info_list(parSampleFile6, meta, "decontX_contamination", "decontX", is_character=FALSE)
}

meta$DBT<-"singlet"
if(file.exists(parSampleFile3)){
  meta = fill_meta_info_list(parSampleFile3, meta, "doubletFinder_doublet_label_resolution_1.5", "DF")
  validation_columns<-c(validation_columns, "DF")

  meta = fill_meta_info_list(parSampleFile3, meta, c("scDblFinder_doublet_call", "scDblFinder_class"), "SDF")
  validation_columns<-c(validation_columns, "SDF")

  meta = fill_meta_info_list(parSampleFile3, meta, "scds_hybrid_call", "scds")
  if(is.logical(meta$scds)){
    meta$scds = ifelse(meta$scds, "Doublet", "Singlet")
  }
  validation_columns<-c(validation_columns, "scds")

  if(!has_decontX){
    meta = fill_meta_info_list(parSampleFile3, meta, "decontX_contamination", "decontX", is_character=FALSE)    
    has_decontX = TRUE
  }
}

if(exists("parSampleFile4")){
  meta = fill_meta_info_list(parSampleFile4, meta, "signacx_CellStates", "SignacX")
  validation_columns<-c(validation_columns, "SignacX")
}

if(exists('parSampleFile5')){
  meta = fill_meta_info_list(parSampleFile5, meta, "SingleR_labels", "SingleR")
  validation_columns<-c(validation_columns, "SingleR")
}

if(exists('parSampleFile7')){
  meta = fill_meta_info_list(parSampleFile7, meta, "Azimuth_finest", "Azimuth")
  validation_columns<-c(validation_columns, "Azimuth")
}

saveRDS(meta, paste0(file_prefix, ".meta.rds"))

obj<-read_object(parFile1)
obj@meta.data = meta

if(length(unique(meta$orig.ident)) > 1){
  validation_columns<-c("orig.ident", validation_columns)
}

do_validation<-function(obj, ct_meta, validation_column, file_prefix, pct){
  tbl = as.data.frame(table(ct_meta[,validation_column]))
  tbl=tbl[order(tbl$Freq, decreasing=TRUE),]
  write.csv(tbl, paste0(file_prefix, ".", pct, ".", validation_column, ".csv"), row.names=FALSE)

  ct_obj = get_filtered_obj(obj, validation_column, ct_meta)
  cur_ct=unique(ct_obj@meta.data[,validation_column])

  g<-get_dim_plot_labelby(ct_obj, reduction="umap", label.by=validation_column,  title=paste0(validation_column, " in old UMAP")) + guides(fill=guide_legend(ncol =1))
  ggsave(paste0(file_prefix, ".", pct, ".", validation_column, ".png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")

  validation_cell=paste0(validation_column, "_Cell")

  ct_obj@meta.data = add_column_count(ct_obj@meta.data, validation_column, validation_cell)
  g<-get_bubble_plot( ct_obj, 
                      assay=cur_assay,
                      group.by=validation_cell, 
                      bubblemap_file = bubblemap_file, 
                      species=species)
  ggsave(paste0(file_prefix, ".", pct, ".", validation_column, ".bubble.png"), g, width=get_dot_width(g), height=get_dot_height(ct_obj, validation_column), units="px", dpi=300, bg="white")

  return(cur_ct)
}

cts = unique(meta[,celltype_column])

ct = cts[1]
for(ct in cts){
  pct = celltype_to_filename(ct)
  ct_meta = meta[meta[,celltype_column] == ct,]

  bar_file=paste0(file_prefix, ".", pct, ".png")
  g<-get_barplot( ct_meta=ct_meta, 
                  bar_file=bar_file,
                  cluster_name=celltype_cluster_column, 
                  validation_columns=validation_columns,
                  calc_height_per_cluster=200, 
                  calc_width_per_cell=50)

  if(has_decontX){
    cells = rownames(ct_meta)
    ct_obj = subset(obj, cells=cells)

    g1<-MyFeaturePlot(ct_obj, features = "decontX") + xlab("") + theme_bw3(TRUE) + theme(aspect.ratio=1) + ggtitle("")
    g2<-VlnPlot(ct_obj, features = "decontX", group.by=celltype_cluster_column) + xlab("") + theme_bw3(TRUE) + ggtitle("") + NoLegend()
    if("DF" %in% colnames(ct_obj@meta.data)){
      g2$data$DF = ct_obj@meta.data[rownames(g2$data), "DF"]
      g2<-g2 + facet_grid(DF~.)
    }
    g<-g1+g2+plot_layout(design="ABBB")
    ggsave(paste0(file_prefix, ".", pct, ".decontX.png"), g, width=4400, height=1600, units="px", dpi=300, bg="white")
  }

  if("SingleR" %in% validation_columns){
    singleR_ct = do_validation(obj, ct_meta, "SingleR", file_prefix, pct)
  }

  if("SignacX" %in% validation_columns){
    signacX_ct = do_validation(obj, ct_meta, "SignacX", file_prefix, pct)
  }

  if("Azimuth" %in% validation_columns){
    azimuth_ct = do_validation(obj, ct_meta, "Azimuth", file_prefix, pct)
  }

  if(all(c("SignacX", "SingleR") %in% validation_columns)){
    tbl = as.data.frame.matrix(table(ct_meta$SingleR, ct_meta$SignacX))
    tbl=tbl[rownames(tbl) %in% singleR_ct, colnames(tbl) %in% signacX_ct]
    write.csv(tbl, paste0(file_prefix, ".", pct, ".SingleR_SignacX.csv"))
  }
}

writeLines(validation_columns, "validation_columns.txt")
if(file.exists(parFile3)){
  writeLines(parFile3, "iter_png.txt")
}

