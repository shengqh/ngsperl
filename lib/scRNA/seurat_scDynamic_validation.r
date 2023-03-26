rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.final.rds'
parFile2='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose/result/AK6383.meta.rds'
parFile3='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_SignacX/result/AK6383.meta.rds'
parFile4='/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_SingleR/result/AK6383.meta.rds'


setwd('/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires_03_choose_validation/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)

options(future.globals.maxSize= 10779361280)

option_tbl=read.table(parSampleFile1, sep="\t")
myoptions = split(option_tbl$V1, option_tbl$V2)
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

meta$DBT<-"singlet"
if(file.exists(parSampleFile2)){
  sctk_files<-read.table(parSampleFile2, sep="\t")
  sctk_map<-unlist(split(sctk_files$V1, sctk_files$V2))

  cur_name='AP_1'
  for(cur_name in names(sctk_map)){
    sctk_meta_file=sctk_map[cur_name]
    sctk_meta=readRDS(sctk_meta_file)

    cur_names<-paste0(cur_name, "_", rownames(sctk_meta))
    cur_cells = intersect(cur_names, rownames(meta))

    if(length(cur_cells) == 0){
      if("project" %in% colnames(meta)){
        obj_meta = meta[meta$project == cur_name,]
      }else{
        obj_meta = meta[meta$sample == cur_name,]
      }
      obj_meta$original_cells<-gsub(".+_", "", rownames(obj_meta))
      cur_cells = intersect(rownames(sctk_meta), obj_meta$original_cells)
      if(length(cur_cells) == 0){
        stop(paste0("I don't know how to map sctk meta cell ", rownames(sctk_meta)[1], " with object meta ", rownames(meta)[1]))
      }else{
        cell_map = unlist(split(obj_meta$original_cells, rownames(obj_meta)))
        meta[names(cell_map), "DBT"] = tolower(as.character(sctk_meta[cell_map,doublet_column]))
      }
    }else{
      rownames(sctk_meta) = cur_names
      meta[cur_cells, "DBT"] = tolower(as.character(sctk_meta[cur_cells,doublet_column]))
    }
  }
  validation_columns<-c(validation_columns, "DBT")
}

saveRDS(meta, paste0(outFile, ".meta.rds"))

draw_figure<-function(outFile, meta, celltype_column, celltype_cluster_column, validation_columns){
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
  }
}

validation_columns<-c("orig.ident", validation_columns)
draw_figure(outFile, meta, celltype_column, celltype_cluster_column, validation_columns)

writeLines(validation_columns, "validation_columns.txt")
if(exists("parFile5")){
  writeLines(parFile5, "iter_png.txt")
}

