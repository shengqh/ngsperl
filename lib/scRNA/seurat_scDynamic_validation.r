rm(list=ls()) 
outFile='P5798'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parSampleFile6='fileList6.txt'
parFile1='/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/seurat_sct2_merge/result/P5798.final.rds'
parFile2='/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/seurat_sct2_merge_dr0.1_1_call/result/P5798.scDynamic.meta.rds'
parFile3='/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/seurat_sct2_merge_dr0.1_1_call/result/P5798.iter_png.csv'


setwd('/nobackup/brown_lab/projects/20231114_scRNA_5798_human_liver_redo/seurat_sct2_merge_dr0.1_1_call_validation/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)

options(future.globals.maxSize= 10779361280)

myoptions = read_file_map(parSampleFile1, do_unlist = FALSE)
doublet_column = myoptions$doublet_column
celltype_column = myoptions$celltype_column

file_dir=paste0(outFile, ".dynamic_call_validation")
dir.create(file_dir, showWarnings=FALSE)
file_prefix=file.path(file_dir, outFile)

meta<-readRDS(parFile2)

if(celltype_column == "seurat_cell_type"){
  celltype_cluster_column = "seurat_clusters"
}else if (paste0("seurat_", celltype_column) %in% colnames(meta)){
  celltype_cluster_column = "seurat_clusters"
}else{
  celltype_cluster_column = paste0(celltype_column, "_clusters")
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

saveRDS(meta, paste0(file_prefix, ".meta.rds"))

obj<-read_object(parFile1)
obj@meta.data = meta

get_filtered_obj<-function(obj, ct_meta, filter_column){
  ct_tbl = table(ct_meta[,filter_column])
  ct_tbl = ct_tbl / sum(ct_tbl)
  ct_tbl = ct_tbl[ct_tbl > 0.01]
  cur_meta = ct_meta[ct_meta[,filter_column] %in% names(ct_tbl),]
  cells = rownames(cur_meta)
  ct_obj = subset(obj, cells=cells)
  return(ct_obj)
}

draw_figure<-function(file_prefix, meta, celltype_column, celltype_cluster_column, validation_columns, has_decontX, obj){
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
      ct_obj = get_filtered_obj(obj, ct_meta, "SingleR")

      g<-get_dim_plot_labelby(ct_obj, reduction="umap", label.by="SingleR",  title="SingleR in old UMAP") + guides(fill=guide_legend(ncol =1))
      ggsave(paste0(file_prefix, ".", pct, ".SingleR.png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")
    }

    if("SignacX" %in% validation_columns){
      ct_obj = get_filtered_obj(obj, ct_meta, "SignacX")

      g<-get_dim_plot_labelby(ct_obj, reduction="umap", label.by="SignacX",  title="SignacX in old UMAP") + guides(fill=guide_legend(ncol =1))
      ggsave(paste0(file_prefix, ".", pct, ".SignacX.png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")
    }
  }
}

if(length(unique(meta$orig.ident)) > 1){
  validation_columns<-c("orig.ident", validation_columns)
}

draw_figure(file_prefix, meta, celltype_column, celltype_cluster_column, validation_columns, has_decontX, obj)

writeLines(validation_columns, "validation_columns.txt")
if(file.exists(parFile3)){
  writeLines(parFile3, "iter_png.txt")
}
