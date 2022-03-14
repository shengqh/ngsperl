
clonos<-read.csv(parFile1, stringsAsFactors=F)
clonos<-clonos[order(clonos$frequency, decreasing=T),]
clono_cells<-unique(clonos[,c("clonotype_id", "cells")])

ct<-read.csv(parFile2, row.names=1)
ct$ident_celltype<-paste0(ct$orig.ident, ":", ct$seurat_cluster, ":", ct$cellactivity_clusters)

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

library(reshape2)
clon_ic_mdf<-data.frame(acast(clon_ic_df, clonotype_id~ctv, value.var="Freq", fill=0), check.names=F)
clon_ic_mdf$clonotype_id<-rownames(clon_ic_mdf)

final <-
  merge(clonos, clon_ic_mdf,
        by = "clonotype_id",
        all.x = T,
        sort = F)

write.csv(
  final,
  paste0(outFile, ".clonotype_cluster.csv"),
  row.names = F,
  quote = F
)
