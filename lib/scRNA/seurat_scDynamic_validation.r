rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parSampleFile6='fileList6.txt'
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_sct_merge/result/combined.final.rds'
parFile2='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_sct_merge_dr0.5_01_call/result/combined.scDynamic.meta.rds'
parFile3='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_sct_merge_dr0.5_01_call/result/combined.iter_png.csv'


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_sct_merge_dr0.5_01_call_validation/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)

options(future.globals.maxSize= 10779361280)

myoptions = read_file_map(parSampleFile1, do_unlist = FALSE)
doublet_column = myoptions$doublet_column
celltype_column = myoptions$celltype_column

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
  meta = fill_meta_info_list(parSampleFile3, meta, "scDblFinder_class", "SDF")
  validation_columns<-c(validation_columns, "SDF")
  meta = fill_meta_info_list(parSampleFile3, meta, "scds_hybrid_call", "scds")
  meta$scds = ifelse(meta$scds, "Doublet", "Singlet")
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

saveRDS(meta, paste0(outFile, ".meta.rds"))

if(has_decontX){
  obj<-read_object(parFile1)
  obj@meta.data = meta
}else{
  obj<-NA
}

draw_figure<-function(outFile, meta, celltype_column, celltype_cluster_column, validation_columns, has_decontX, obj){
  cts = unique(meta[,celltype_column])

  ct = cts[1]
  for(ct in cts){
    pct = celltype_to_filename(ct)
    ct_meta = meta[meta[,celltype_column] == ct,]

    alltbl=NULL

    col_name="SignacX"
    for(col_name in validation_columns){
      tbl = data.frame(table(ct_meta[,celltype_cluster_column], ct_meta[,col_name]))
      v1 = as.numeric(as.character(tbl$Var1))
      if(all(is.na(v1))){
        v1 = as.character(tbl$Var1)
      }
      tbl$Var1 = v1
      tbl$Category=col_name

      alltbl<-rbind(alltbl, tbl)
    }

    g<-ggplot(alltbl, aes(Var2, Freq, fill=Var2)) + geom_bar(width=0.5, stat = "identity") + facet_grid(Var1~Category, scales = "free", space='free_x') + theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend()

    height = max(500, length(unique(alltbl$Var1)) * 200) + 500
    width = max(1000, length(unique(alltbl$Var2)) * 50) + 400
    png(paste0(outFile, ".", pct, ".png"), width=width, height=height, res=300)
    print(g)
    dev.off()

    if(has_decontX){
      g1<-MyFeaturePlot(obj, features = "decontX") + xlab("") + theme_bw3(TRUE) + theme(aspect.ratio=1) + ggtitle("")
      g2<-VlnPlot(obj, features = "decontX", group.by=celltype_cluster_column) + xlab("") + theme_bw3(TRUE) + ggtitle("") + NoLegend()
      if("DF" %in% colnames(obj@meta.data)){
        g2$data$DF = obj@meta.data[rownames(g2$data), "DF"]
        g2<-g2 + facet_grid(DF~.)
      }
      g<-g1+g2+plot_layout(design="ABBB")
      png(paste0(outFile, ".", pct, ".decontX.png"), width=4400, height=1600, res=300)
      print(g)
      dev.off()
    }
  }
}

if(length(unique(meta$orig.ident)) > 1){
  validation_columns<-c("orig.ident", validation_columns)
}

draw_figure(outFile, meta, celltype_column, celltype_cluster_column, validation_columns, has_decontX, obj)

writeLines(validation_columns, "validation_columns.txt")
if(file.exists(parFile3)){
  writeLines(parFile3, "iter_png.txt")
}
