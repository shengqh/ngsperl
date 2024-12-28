rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parSampleFile6='fileList6.txt'
parSampleFile8='fileList8.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20241224_combined_scRNA_hg38_cellbender_fastmnn/cellbender_raw_qc_report/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(reshape2)
library(SingleCellExperiment)

options(future.globals.maxSize= 10779361280)

myoptions = read_file_map(parSampleFile1, do_unlist=FALSE)
doublet_column = myoptions$doublet_column
celltype_cluster_column = myoptions$celltype_column
celltype_column="orig.ident"

rmd_ext=gsub(".html","",myoptions$rmd_ext)
detail_folder=paste0(myoptions$prefix, rmd_ext)
dir.create(detail_folder, showWarnings = FALSE)
detail_prefix=file.path(detail_folder, myoptions$prefix)

obj_map<-read_file_map(parSampleFile2)

sample_names=names(obj_map)

validation_columns=c()

has_sctk<-file.exists(parSampleFile3)
if(has_sctk){
  sctk_map<-read_file_map(parSampleFile3)
  validation_columns<-c(validation_columns, "sctk_DF", "sctk_SDF", "sctk_scds")
}

has_signacx<-exists('parSampleFile4')
if(has_signacx){
  signacx_map<-read_file_map(parSampleFile4)
  validation_columns<-c(validation_columns, "SignacX")
}

has_singler<-exists('parSampleFile5')
if(has_singler){
  singler_map<-read_file_map(parSampleFile5)
  validation_columns<-c(validation_columns, "SingleR")
}

has_decontX<-exists('parSampleFile6')
if(has_decontX){
  decontX_map<-read_file_map(parSampleFile6)
  validation_columns<-c(validation_columns, "decontX > 0.25")
}

has_doublet_finder<-exists('parSampleFile7')
if(has_doublet_finder){
  doublet_finder_map<-read_file_map(parSampleFile7)
}

has_azimuth<-exists('parSampleFile8')
if(has_azimuth){
  azimuth_map<-read_file_map(parSampleFile8)
  validation_columns<-c(validation_columns, "Azimuth")
}

draw_figure<-function(sample_prefix, cur_meta, cur_validation_columns){
  alltbl=NULL

  col_name="SignacX"
  for(col_name in cur_validation_columns){
    tbl = data.frame(table(cur_meta[,"seurat_cell_type"], cur_meta[,col_name]))
    tbl$Category=col_name
    alltbl<-rbind(alltbl, tbl)
  }

  levels(alltbl$Var1)<-levels(cur_meta$seurat_cell_type)

  g<-ggplot(alltbl, aes(Var2, Freq, fill=Var2)) + 
    geom_bar(width=0.5, stat = "identity") + 
    facet_grid(Var1~Category, scales = "free", space='free_x') + 
    theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend() + 
    theme(strip.text.y.right = element_text(angle = 0, hjust = 0),
          strip.text.x.top = element_text(angle = 90, hjust = 0))

  height = max(800, length(unique(alltbl$Var1)) * 150) + 500
  width = max(1000, length(unique(alltbl$Var2)) * 50) + 1000

  png(paste0(sample_prefix, ".validation.png"), width=width, height=height, res=300)
  print(g)
  dev.off()
}

meta<-NULL
stats_df<-NULL
ct_tb<-NULL
sample_name = sample_names[1]
for(sample_name in sample_names){
  cat("read", sample_name, "...\n")
  obj_file = obj_map[[sample_name]]

  sample_prefix = file.path(detail_folder, sample_name)

  object.list<-readRDS(obj_file)
  stats_df<-rbind(stats_df, object.list[[sample_name]]$preprocess)

  obj = object.list[[sample_name]]$obj
  cur_meta<-obj@meta.data

  ct_df<-data.frame(table(cur_meta$cell_type))
  ct_df$Sample<-sample_name

  ct_tb<-rbind(ct_tb, ct_df)

  cur_validation_columns = validation_columns;
  
  df_column = ""
  
  if(has_sctk){
    sctk_file = sctk_map[[sample_name]]
    sctk_meta = readRDS(sctk_file)

    cur_meta = fill_meta_info(sample_name, sctk_meta, cur_meta, "doubletFinder_doublet_label_resolution_1.5", "sctk_DF")
    df_column = "sctk_DF"

    cur_meta = fill_meta_info(sample_name, sctk_meta, cur_meta, c("scDblFinder_doublet_call", "scDblFinder_class"), "sctk_SDF")

    cur_meta = fill_meta_info(sample_name, sctk_meta, cur_meta, "scds_hybrid_call", "sctk_scds")

    if(is.logical(cur_meta$scds)){
      cur_meta$scds = ifelse(cur_meta$scds, "Doublet", "Singlet")
    }
  }
  
  if(has_signacx){
    signacx_file = signacx_map[[sample_name]]
    signacx_meta = readRDS(signacx_file)
    cur_meta = fill_meta_info(sample_name, signacx_meta, cur_meta, "signacx_CellStates", "SignacX")
  }
  
  if(has_singler){
    singler_file = singler_map[[sample_name]]
    singler_meta = readRDS(singler_file)
    cur_meta = fill_meta_info(sample_name, singler_meta, cur_meta, "SingleR_labels", "SingleR")
  }
  
  if(has_azimuth){
    azimuth_file = azimuth_map[[sample_name]]
    azimuth_meta = readRDS(azimuth_file)
    cur_meta = fill_meta_info(sample_name, azimuth_meta, cur_meta, "Azimuth_finest", "Azimuth")
  }
  
  if(has_doublet_finder){
    df_file = doublet_finder_map[[sample_name]]
    df_meta = readRDS(df_file)

    df_option_file = gsub(".meta.rds", ".options.csv", df_file)
    df_option = read.csv(df_option_file)
    
    idx=1
    for(idx in c(1:nrow(df_option))){
      df_rate = df_option$doublet_rate[idx]
      df_label = df_option$label[idx]
      df_column = paste0("DF_", round(df_rate, 3))
      cur_meta = fill_meta_info(sample_name, df_meta, cur_meta, df_label, df_column)
      cur_validation_columns = c(cur_validation_columns, df_column)
    }
  }

  if(has_decontX){
    decontX_file = decontX_map[[sample_name]]
    decontX_meta = readRDS(decontX_file)
    cur_meta = fill_meta_info(sample_name, decontX_meta, cur_meta, "decontX_contamination", "decontX_contamination", is_character = FALSE)
    cur_meta[,"decontX > 0.25"] = cur_meta[,"decontX_contamination"] > 0.25
    obj@meta.data = cur_meta
    
    g1<-MyFeaturePlot(obj, features = "decontX_contamination") + xlab("") + theme_bw3(TRUE) + theme(aspect.ratio=1) + ggtitle("")
    g2<-VlnPlot(obj, features = "decontX_contamination", group.by="seurat_cell_type") + xlab("") + theme_bw3(TRUE)  + ggtitle("") + NoLegend()
    if(df_column != ""){
      g2$data$DF = obj@meta.data[rownames(g2$data), df_column]
      g2<-g2 + facet_grid(DF~.)
    }
    g<-g1+g2+plot_layout(design="ABBB")
    ggsave(paste0(sample_prefix, ".decontX.png"), width=4400, height=1600, dpi=300, units="px", bg="white")
  }
  
  rownames(cur_meta)<-paste0(sample_name, "_", rownames(cur_meta))
  cat("save meta ...\n")
  saveRDS(cur_meta, paste0(sample_prefix, ".meta.rds"))

  cat("draw_figure ...\n")
  draw_figure(sample_prefix, cur_meta, cur_validation_columns)
}

stats_df<-stats_df[,colnames(stats_df) != "tringsAsFactors"]
if(all(colnames(stats_df)[1:2] == "sample")){
  stats_df=stats_df[,c(2:ncol(stats_df))]
}
write.csv(stats_df, file.path(detail_folder, "sample_summary.csv"), row.names=F)

ct_tb<-acast(ct_tb, "Sample~Var1",  value.var="Freq", fill=0)
write.csv(ct_tb, file.path(detail_folder, "sample_celltype.csv"), row.names=TRUE)
writeLines(validation_columns, file.path(detail_folder,"validation_columns.txt"))
