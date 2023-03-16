rm(list=ls()) 
outFile='crs'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/nobackup/h_turner_lab/shengq2/20230314_7114_8822_scRNA_hg38/seurat_sct_merge/result/crs.final.rds'
parFile2='/nobackup/h_turner_lab/shengq2/20230314_7114_8822_scRNA_hg38/seurat_sct_merge_dr0.5_01_call/result/crs.scDynamic.meta.rds'
parFile3='/nobackup/h_turner_lab/shengq2/20230314_7114_8822_scRNA_hg38/seurat_sct_merge_SignacX/result/crs.meta.rds'
parFile4='/nobackup/h_turner_lab/shengq2/20230314_7114_8822_scRNA_hg38/seurat_sct_merge_singleR/result/crs.meta.rds'


setwd('/nobackup/h_turner_lab/shengq2/20230314_7114_8822_scRNA_hg38/seurat_sct_merge_dr0.5_01_call_validation/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)

options(future.globals.maxSize= 10779361280)

option_tbl=read.table(parSampleFile1, sep="\t")
myoptions = split(option_tbl$V1, option_tbl$V2)
doublet_column = myoptions$doublet_column

meta<-readRDS(parFile2)

validation_columns=c()

if(file.exists(parFile3)){
  singacx=readRDS(parFile3)
  meta$SignacX=singacx[rownames(meta),"signacx_CellStates"]
  validation_columns<-c(validation_columns, "SignacX")
}

if(exists('parFile4')){
  if(file.exists(parFile4)){
    singleR=readRDS(parFile4)
    meta$SingleR=singleR[rownames(meta),"SingleR_labels"]
    validation_columns<-c(validation_columns, "SingleR")
  }
}

meta$DoubletStatus<-"singlet"
if(file.exists(parSampleFile2)){
  sctk_files<-read.table(parSampleFile2, sep="\t")
  sctk_map<-unlist(split(sctk_files$V1, sctk_files$V2))

  cur_name='T6_KS_030822_RC_E'
  for(cur_name in names(sctk_map)){
    cur_meta_file=sctk_map[cur_name]
    cur_meta=readRDS(cur_meta_file)

    rownames(cur_meta)<-paste0(cur_name, "_", rownames(cur_meta))
    cur_cells = intersect(rownames(cur_meta), rownames(meta))
    meta[cur_cells, "DoubletStatus"] = tolower(as.character(cur_meta[cur_cells,doublet_column]))
  }
  validation_columns<-c(validation_columns, "DoubletStatus")
}

saveRDS(meta, paste0(outFile, ".meta.rds"))

col_name="SignacX"
draw_figure<-function(outFile, meta, col_name){
  cts = unique(meta$layer4)
  ct = cts[1]
  for(ct in cts){
    pct = celltype_to_filename(ct)
    ct_meta = subset(meta, layer4 == ct)
    tbl = data.frame(table(ct_meta$layer4_clusters, ct_meta[,col_name]))
    tbl$Var1 = as.numeric(as.character(tbl$Var1))

    g<-ggplot(tbl, aes(Var2, Freq, fill=Var2)) + geom_bar(width=0.5, stat = "identity") + facet_grid(Var1~., scales = "free_y") + theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend()

    height = max(1000, length(unique(tbl$Var1)) * 200) + 500
    width = max(1000, length(unique(tbl$Var2)) * 100) + 400
    png(paste0(outFile, ".", pct, ".", col_name, ".png"), width=width, height=height, res=300)
    print(g)
    dev.off()
  }
}

for(col_name in validation_columns){
  draw_figure(outFile, meta, col_name)
}
