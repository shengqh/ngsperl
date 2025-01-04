rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parSampleFile7='fileList7.txt'
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.final.rds'
parFile2='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.meta.rds'
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose_validation/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)

options(future.globals.maxSize= 10779361280)

myoptions = read_file_map(parSampleFile1, do_unlist = FALSE)
reduction = myoptions$reduction
doublet_column = myoptions$doublet_column
celltype_column = myoptions$celltype_column
bubblemap_file = myoptions$bubblemap_file
cur_assay = ifelse(is_one(myoptions$by_sctransform), "SCT", "RNA")
species = myoptions$species
create_clusters = is_one(myoptions$create_clusters)
summary_layer=myoptions$summary_layer

cat("reduction: ", reduction, "\n")
cat("doublet_column: ", doublet_column, "\n")
cat("celltype_column: ", celltype_column, "\n")

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
    meta[,"decontX > 0.25"] = meta[,"decontX"] > 0.25
    has_decontX = TRUE
    validation_columns<-c(validation_columns, "decontX > 0.25")
  }
}

if(!("SignacX" %in% colnames(meta))){
  if(exists("parSampleFile4")){
    meta = fill_meta_info_list(parSampleFile4, meta, "signacx_CellStates", "SignacX")
    validation_columns<-c(validation_columns, "SignacX")
  }
}else{
  validation_columns<-c(validation_columns, "SignacX")
}

if(!("SingleR" %in% colnames(meta))){
  if(exists('parSampleFile5')){
    meta = fill_meta_info_list(parSampleFile5, meta, "SingleR_labels", "SingleR")
    validation_columns<-c(validation_columns, "SingleR")
  }
}else{
  validation_columns<-c(validation_columns, "SingleR")
}

if(!("Azimuth" %in% colnames(meta))) {
  if(exists('parSampleFile7')){
    meta = fill_meta_info_list(parSampleFile7, meta, "Azimuth_finest", "Azimuth")
    validation_columns<-c(validation_columns, "Azimuth")
  }
}else{
  validation_columns<-c(validation_columns, "Azimuth")
}
saveRDS(meta, paste0(file_prefix, ".meta.rds"))

cat("reading object\n")
obj<-read_object(parFile1)
obj@meta.data = meta

if(length(unique(meta$orig.ident)) > 1){
  validation_columns<-c("orig.ident", validation_columns)
}

do_validation<-function(obj, ct_meta, validation_column, file_prefix, pct, reduction="umap", summary_layer="layer4"){
  tbl = as.data.frame(table(ct_meta[,validation_column]))
  tbl=tbl[order(tbl$Freq, decreasing=TRUE),]
  write.csv(tbl, paste0(file_prefix, ".", pct, ".", validation_column, ".csv"), row.names=FALSE)

  ct_obj = get_filtered_obj(obj, validation_column, ct_meta)
  cur_ct=unique(ct_obj@meta.data[,validation_column])

  g<-get_dim_plot_labelby(ct_obj, reduction=reduction, label.by=validation_column,  title=paste0(validation_column, " in UMAP")) + guides(fill=guide_legend(ncol =1))
  ggsave(paste0(file_prefix, ".", pct, ".", validation_column, ".png"), g, width=2000, height=1200, units="px", dpi=300, bg="white")

  if(summary_layer %in% colnames(obj@meta.data)){
    g<-get_sub_bubble_plot(
      obj=obj, 
      obj_res=summary_layer, 
      subobj=ct_obj, 
      subobj_res=validation_column, 
      bubblemap_file=bubblemap_file, 
      add_num_cell=TRUE, 
      species=NULL, 
      assay="RNA")
  }else{
    validation_cell=paste0(validation_column, "_Cell")
    
    ct_obj@meta.data = add_column_count(ct_obj@meta.data, validation_column, validation_cell)
    g<-get_bubble_plot( ct_obj, 
                        assay=cur_assay,
                        group.by=validation_cell, 
                        bubblemap_file = bubblemap_file, 
                        species=species)
  }

  ggsave(paste0(file_prefix, ".", pct, ".", validation_column, ".bubble.png"), g, width=get_dot_width(g), height=get_dot_height(ct_obj, validation_column), units="px", dpi=300, bg="white")

  return(cur_ct)
}

if(!is.null(summary_layer)){
  scts = unique(meta[,summary_layer])
  sct=scts[1]
  for(sct in scts){
    cat("dimplot of", sct, "\n")
    pct=celltype_to_filename(sct)
    sct_obj=subset(obj, cells=rownames(meta)[meta[,summary_layer]==sct])
    g=get_dim_plot( sct_obj, 
                    reduction=reduction, 
                    group.by=celltype_cluster_column,
                    label.by=celltype_column, 
                    title=sct) + guides(fill=guide_legend(ncol =1))
    ggsave( paste0(file_prefix, ".", pct, ".umap.png"), 
            g, 
            width=2000, 
            height=1200, 
            units="px", 
            dpi=300, 
            bg="white")
  }
}

cts = levels(meta[,celltype_column])
ct = cts[1]
for(ct in cts){
  cat(ct, "\n")
  pct = celltype_to_filename(ct)
  ct_meta = meta[meta[,celltype_column] == ct,]

  cat("  barplot\n")
  bar_file=paste0(file_prefix, ".", pct, ".png")
  g<-get_barplot( ct_meta=ct_meta, 
                  bar_file=bar_file,
                  cluster_name=celltype_cluster_column, 
                  validation_columns=validation_columns,
                  calc_height_per_cluster=200, 
                  calc_width_per_cell=50)

  cat("  highlight plot\n")
  g = Meta_Highlight_Plot(seurat_object = obj, 
                          meta_data_column = celltype_column,
                          meta_data_highlight = ct, 
                          highlight_color = c("firebrick"),
                          background_color = "lightgray",
                          reduction="umap") + 
                          theme(legend.position="top")
  if(celltype_column == "layer4"){
    ggsave(paste0(file_prefix, ".", pct, ".umap.png"), g, width=1200, height=1300, units="px", dpi=300, bg="white")
  }else{
    layer4_ct=names(table(as.character(ct_meta$layer4)))[1]
    layer4_cells=rownames(meta)[meta$layer4==layer4_ct]
    layer4_obj=subset(obj, cells=layer4_cells)
    g2 = Meta_Highlight_Plot(seurat_object = layer4_obj, 
                            meta_data_column = celltype_column,
                            meta_data_highlight = ct, 
                            highlight_color = c("firebrick"),
                            background_color = "lightgray",
                            reduction="subumap") + 
                            theme(legend.position="top")
    g=g+g2+plot_layout(design="AB")
    ggsave(paste0(file_prefix, ".", pct, ".umap.png"), g, width=2400, height=1300, units="px", dpi=300, bg="white")
  }

  if(has_decontX){
    cells = rownames(ct_meta)
    ct_obj = subset(obj, cells=cells)

    cat("  decontX plot\n")
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
    cat("  SingleR plot\n")
    singleR_ct = do_validation( obj, 
                                ct_meta, 
                                "SingleR", 
                                file_prefix, 
                                pct,
                                reduction=reduction,
                                summary_layer=summary_layer)
  }

  if("SignacX" %in% validation_columns){
    cat("  SignacX plot\n")
    signacX_ct = do_validation( obj, 
                                ct_meta, 
                                validation_column = "SignacX", 
                                file_prefix, 
                                pct,
                                reduction=reduction,
                                summary_layer=summary_layer)
  }

  if("Azimuth" %in% validation_columns){
    cat("  Azimuth plot\n")
    azimuth_ct = do_validation( obj, 
                                ct_meta, 
                                "Azimuth", 
                                file_prefix, 
                                pct,
                                reduction=reduction,
                                summary_layer=summary_layer)
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
