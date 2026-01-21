rm(list=ls()) 
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/20260120_silhouetee_pairwise/silhouette_pairwise_merge/result')

### Parameter setting end ###

library(Seurat)

file_tbl=read.table(parSampleFile1, header=FALSE, sep="\t")

myoptions=read.table(parSampleFile2, header=FALSE, sep="\t")

obj_file=myoptions[myoptions$V2==outFile, 1]
cluster=myoptions[myoptions$V2=="cluster", 1]

cat("Reading object file:", obj_file, "\n")
obj=readRDS(obj_file)

sil.data=data.frame(cluster=FetchData(obj, vars = cluster)[, 1], neighbor=100, sil_width=100)
rownames(sil.data)=colnames(obj)

cur_file=file_tbl$V1[2]
for(cur_file in file_tbl$V1){
  if(!file.exists(cur_file)) {
    cat("File does not exist:", cur_file, "\n")
    next
  }

  cat("Reading silhouette file:", cur_file, "\n")
  sil_list=readRDS(cur_file)
  i = sil_list$i
  j = sil_list$j
  sil = sil_list$silhouette

  #combine sil with sil.data only if the sil_width in sil less than sil_width in sil.data

  stopifnot(all(rownames(sil) %in% rownames(sil.data)))
 
  cur_sil = sil.data[rownames(sil),]
 
  update_cells = rownames(sil)[sil$sil_width < cur_sil$sil_width]
  if(length(update_cells) > 0){
    cat("Updating silhouette width for", length(update_cells), "cells\n")
    sil.data[update_cells, "sil_width"] = sil[update_cells, "sil_width"]
    sil.data[update_cells, "neighbor"] = sil[update_cells, "neighbor"]
  }
}

sil.data$closest <- factor(ifelse(sil.data$sil_width > 0, sil.data$cluster, sil.data$neighbor))
sil.data_rds <- paste0(outFile, ".silhouetee.rds")
saveRDS(sil.data, sil.data_rds)
