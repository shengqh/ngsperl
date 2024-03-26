rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parSampleFile7='fileList7.txt'
parFile1='/data/h_gelbard_lab/projects/20240320_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.2_3_choose/result/iSGS_cell_atlas.final.rds'
parFile2='/data/h_gelbard_lab/projects/20240320_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.2_3_choose/result/iSGS_cell_atlas.meta.rds'
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240320_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.2_3_choose_validation/result')

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

file_dir=paste0(outFile, gsub(".html", "", myoptions$rmd_ext))
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
    ct_obj = get_filtered_obj(obj, "SingleR", ct_meta)
    singleR_ct=unique(ct_obj$SingleR)

    g<-get_dim_plot_labelby(ct_obj, reduction="umap", label.by="SingleR",  title="SingleR in old UMAP") + guides(fill=guide_legend(ncol =1))
    ggsave(paste0(file_prefix, ".", pct, ".SingleR.png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")

    ct_obj@meta.data = add_column_count(ct_obj@meta.data, "SingleR", "SingleR_Cell")
    g<-get_bubble_plot( ct_obj, 
                        assay=cur_assay,
                        group.by="SingleR_Cell", 
                        bubblemap_file = bubblemap_file, 
                        species=species)
    ggsave(paste0(file_prefix, ".", pct, ".SingleR.bubble.png"), g, width=get_dot_width(g), height=get_dot_height(ct_obj, "SingleR"), units="px", dpi=300, bg="white")
  }

  if("SignacX" %in% validation_columns){
    ct_obj = get_filtered_obj(obj, "SignacX", ct_meta)
    signacX_ct=unique(ct_obj$SignacX)

    g<-get_dim_plot_labelby(ct_obj, reduction="umap", label.by="SignacX",  title="SignacX in old UMAP") + guides(fill=guide_legend(ncol =1))
    ggsave(paste0(file_prefix, ".", pct, ".SignacX.png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")

    ct_obj@meta.data = add_column_count(ct_obj@meta.data, "SignacX", "SignacX_Cell")
    g<-get_bubble_plot( ct_obj, 
                        assay=cur_assay,
                        group.by="SignacX_Cell", 
                        bubblemap_file = bubblemap_file, 
                        species=species)
    ggsave(paste0(file_prefix, ".", pct, ".SignacX.bubble.png"), g, width=get_dot_width(g), height=get_dot_height(ct_obj, "SignacX"), units="px", dpi=300, bg="white")
  }

  if("Azimuth" %in% validation_columns){
    ct_obj = get_filtered_obj(obj, "Azimuth", ct_meta)

    g<-get_dim_plot_labelby(ct_obj, reduction="umap", label.by="Azimuth",  title="Azimuth in old UMAP") + guides(fill=guide_legend(ncol =1))
    ggsave(paste0(file_prefix, ".", pct, ".Azimuth.png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")

    ct_obj@meta.data = add_column_count(ct_obj@meta.data, "Azimuth", "Azimuth_Cell")
    g<-get_bubble_plot( ct_obj, 
                        assay=cur_assay,
                        group.by="Azimuth_Cell", 
                        bubblemap_file = bubblemap_file, 
                        species=species)
    ggsave(paste0(file_prefix, ".", pct, ".Azimuth.bubble.png"), g, width=get_dot_width(g), height=get_dot_height(ct_obj, "Azimuth"), units="px", dpi=300, bg="white")
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
