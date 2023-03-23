rm(list=ls()) 
outFile='P8256'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/scratch/cqs/shengq2/ravi_shah_projects/20230319_validate_code/seurat_merge/result/P8256.final.rds'
parFile2='/scratch/cqs/shengq2/ravi_shah_projects/20230319_validate_code/seurat_merge_dr0.5_01_call/result/P8256.scDynamic.meta.rds'
parFile3='/scratch/cqs/shengq2/ravi_shah_projects/20230319_validate_code/seurat_merge_SignacX/result/P8256.meta.rds'
parFile4='/scratch/cqs/shengq2/ravi_shah_projects/20230319_validate_code/seurat_merge_SingleR/result/P8256.meta.rds'
parFile5='/scratch/cqs/shengq2/ravi_shah_projects/20230319_validate_code/seurat_merge_dr0.5_01_call/result/P8256.iter_png.csv'


setwd('/scratch/cqs/shengq2/ravi_shah_projects/20230319_validate_code/seurat_merge_dr0.5_01_call_validation/result')

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

  cur_name='SB_1'
  for(cur_name in names(sctk_map)){
    sctk_meta_file=sctk_map[cur_name]
    sctk_meta=readRDS(sctk_meta_file)

    cur_names<-paste0(cur_name, "_", rownames(sctk_meta))
    cur_cells = intersect(cur_names, rownames(meta))

    if(length(cur_cells) == 0){
      obj_meta = meta[meta$project == cur_name,]
      obj_meta$original_cells<-gsub(".+_", "", rownames(obj_meta))
      cur_cells = intersect(rownames(sctk_meta), obj_meta$original_cells)
      if(length(cur_cells) == 0){
        stop(paste0("I don't know how to map sctk meta cell ", rownames(sctk_meta)[1], " with object meta ", rownames(meta)[1]))
      }else{
        cell_map = unlist(split(obj_meta$original_cells, rownames(obj_meta)))
        meta[names(cell_map), "Doublet"] = tolower(as.character(sctk_meta[cell_map,doublet_column]))
      }
    }else{
      rownames(sctk_meta) = cur_names
      meta[cur_cells, "Doublet"] = tolower(as.character(sctk_meta[cur_cells,doublet_column]))
    }
  }
  validation_columns<-c(validation_columns, "Doublet")
}

saveRDS(meta, paste0(outFile, ".meta.rds"))

draw_figure<-function(outFile, meta, validation_columns){
  cts = unique(meta$layer4)
  ct = cts[1]
  for(ct in cts){
    pct = celltype_to_filename(ct)
    ct_meta = subset(meta, layer4 == ct)

    alltbl=NULL

    col_name="SignacX"
    for(col_name in validation_columns){
      tbl = data.frame(table(ct_meta$layer4_clusters, ct_meta[,col_name]))
      tbl$Var1 = as.numeric(as.character(tbl$Var1))
      tbl$Category=col_name

      alltbl<-rbind(alltbl, tbl)
    }

    g<-ggplot(alltbl, aes(Var2, Freq, fill=Var2)) + geom_bar(width=0.5, stat = "identity") + facet_grid(Var1~Category, scales = "free", space='free_x') + theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend()

    height = max(1000, length(unique(alltbl$Var1)) * 200) + 500
    width = max(1000, length(unique(alltbl$Var2)) * 50) + 400
    png(paste0(outFile, ".", pct, ".png"), width=width, height=height, res=300)
    print(g)
    dev.off()
  }
}

validation_columns<-c("orig.ident", validation_columns)
draw_figure(outFile, meta, validation_columns)

writeLines(validation_columns, "validation_columns.txt")
writeLines(parFile5, "iter_png.txt")
