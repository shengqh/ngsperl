rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parSampleFile6='fileList6.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230410_combined_scRNA_hg38/raw_qc_sct2_report/result')

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

obj_map<-read_file_map(parSampleFile2)

sample_names=names(obj_map)

validation_columns=c()

has_sctk<-file.exists(parSampleFile3)
if(has_sctk){
  sctk_map<-read_file_map(parSampleFile3)
  validation_columns<-c(validation_columns, "DF", "SDF", "scds")
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
}

draw_figure<-function(sample_name, cur_meta, validation_columns){
  alltbl=NULL

  col_name="SignacX"
  for(col_name in validation_columns){
    tbl = data.frame(table(cur_meta[,"seurat_cell_type"], cur_meta[,col_name]))
    tbl$Category=col_name
    alltbl<-rbind(alltbl, tbl)
  }

  levels(alltbl$Var1)<-levels(cur_meta$seurat_cell_type)

  g<-ggplot(alltbl, aes(Var2, Freq, fill=Var2)) + 
    geom_bar(width=0.5, stat = "identity") + 
    facet_grid(Var1~Category, scales = "free", space='free_x') + 
    theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend() + 
    theme(strip.text.y.right = element_text(angle = 0, hjust = 0))

  height = max(500, length(unique(alltbl$Var1)) * 150) + 500
  width = max(1000, length(unique(alltbl$Var2)) * 50) + 1000

  png(paste0(sample_name, ".validation.png"), width=width, height=height, res=300)
  print(g)
  dev.off()
}

meta<-NULL
stats_df<-NULL
ct_tb<-NULL
sample_name = sample_names[2]
for(sample_name in sample_names){
  cat("read", sample_name, "...\n")
  obj_file = obj_map[[sample_name]]

  object.list<-readRDS(obj_file)
  stats_df<-rbind(stats_df, object.list[[sample_name]]$preprocess)

  obj = object.list[[sample_name]]$obj
  cur_meta<-obj@meta.data

  ct_df<-data.frame(table(cur_meta$cell_type))
  ct_df$Sample<-sample_name

  ct_tb<-rbind(ct_tb, ct_df)

  if(has_sctk){
    sctk_file = sctk_map[[sample_name]]
    sctk_meta = readRDS(sctk_file)

    cur_meta = fill_meta_info(sample_name, sctk_meta, cur_meta, "doubletFinder_doublet_label_resolution_1.5", "DF")
    cur_meta = fill_meta_info(sample_name, sctk_meta, cur_meta, "scDblFinder_class", "SDF")
    cur_meta = fill_meta_info(sample_name, sctk_meta, cur_meta, "scds_hybrid_call", "scds")
    cur_meta$scds = ifelse(cur_meta$scds, "Doublet", "Singlet")
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

  if(has_decontX){
    decontX_file = decontX_map[[sample_name]]
    decontX_meta = readRDS(decontX_file)
    cur_meta = fill_meta_info(sample_name, decontX_meta, cur_meta, "decontX_contamination", "decontX_contamination", is_character = FALSE)
    obj@meta.data = cur_meta

    g1<-MyFeaturePlot(obj, features = "decontX_contamination") + xlab("") + theme_bw3(TRUE) + theme(aspect.ratio=1) + ggtitle("")
    g2<-VlnPlot(obj, features = "decontX_contamination", group.by="seurat_cell_type") + xlab("") + theme_bw3(TRUE)  + ggtitle("") + NoLegend()
    if(has_sctk){
      g2$data$DF = obj@meta.data[rownames(g2$data), "DF"]
      g2<-g2 + facet_grid(DF~.)
    }
    g<-g1+g2+plot_layout(design="ABBB")
    png(paste0(sample_name, ".decontX.png"), width=4400, height=1600, res=300)
    print(g)
    dev.off()
  }

  rownames(cur_meta)<-paste0(sample_name, "_", rownames(cur_meta))
  cat("save meta ...\n")
  saveRDS(cur_meta, paste0(sample_name, ".meta.rds"))

  cat("draw_figure ...\n")
  draw_figure(sample_name, cur_meta, validation_columns)
}

stats_df<-stats_df[,colnames(stats_df) != "tringsAsFactors"]
if(all(colnames(stats_df)[1:2] == "sample")){
  stats_df=stats_df[,c(2:ncol(stats_df))]
}
write.csv(stats_df, "sample_summary.csv", row.names=F)

ct_tb<-acast(ct_tb, "Sample~Var1",  value.var="Freq", fill=0)
write.csv(ct_tb, "sample_celltype.csv", row.names=TRUE)
writeLines(validation_columns, "validation_columns.txt")

