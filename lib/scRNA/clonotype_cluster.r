rm(list=ls()) 
outFile='AG3669'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2='/data/h_gelbard_lab/projects/20220508_scRNA_3669/seurat_merge_03_choose_res/result/AG3669.meta.rds'
parFile3=''
outputDirectory='.'


setwd('/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_04_db_cluster/result')

### Parameter setting end ###

library("tools")
library(reshape2)

file_map = read.table("fileList1.txt", sep="\t")

if(file_ext(parFile2) == "rds"){
  ct<-readRDS(parFile2)
  ct$ident_celltype<-paste0(ct$orig.ident, ":", ct$seurat_cluster, ":", ct$cell_type)
}else{
  ct<-read.csv(parFile2, row.names=1)
  ct$ident_celltype<-paste0(ct$orig.ident, ":", ct$seurat_cluster, ":", ct$cellactivity_clusters)
}

for (parFile1 in file_map$V1){
  samplename <- gsub(".csv", "", basename(parFile1))

  clonos<-read.csv(parFile1, stringsAsFactors=F)
  clonos<-clonos[order(clonos$frequency, decreasing=T),]
  clono_cells<-unique(clonos[,c("clonotype_id", "cells")])

  clon_ic<-apply(clono_cells, 1, function(x){
    clonotype_id=x[['clonotype_id']]
    cells=x[['cells']]
    scells<-unlist(strsplit(cells,';'))
    not_in_cluster<-scells[!(scells %in% rownames(ct))]
    scells<-scells[!(scells %in% not_in_cluster)]
    if(length(scells) > 0){
      ctv<-ct[scells, "ident_celltype"]
      res=data.frame(table(ctv))
      #print(res)
      if(length(not_in_cluster) > 0){
        res<-rbind(res, data.frame("ctv"="not_in_cluster", "Freq"=length(not_in_cluster)))
      }
    }else{
      res<-data.frame("ctv"="not_in_cluster", "Freq"=length(not_in_cluster))
    }
    res$clonotype_id=clonotype_id
    return(res)
  })

  clon_ic_df<-do.call(rbind, clon_ic)
  clon_ic_df$ctv<-factor(clon_ic_df$ctv, levels=sort(unique(as.character(clon_ic_df$ctv))))
  clon_ic_df$clonotype_id<-factor(clon_ic_df$clonotype_id, levels=unique(clon_ic_df$clonotype_id))
  #head(clon_ic_df)

  clon_ic_mdf<-data.frame(acast(clon_ic_df, clonotype_id~ctv, value.var="Freq", fill=0), check.names=F)
  clon_ic_mdf$clonotype_id<-rownames(clon_ic_mdf)

  final <-
    merge(clonos, clon_ic_mdf,
          by = "clonotype_id",
          all.x = T,
          sort = F)

  write.csv(
    final,
    paste0(samplename, ".cluster.csv"),
    row.names = F,
    quote = F
  )
}
