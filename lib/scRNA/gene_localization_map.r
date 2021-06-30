source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(ggpubr)

finalList<-readRDS(parFile1)

obj<-finalList$obj

groups_tbl<-read.table(parSampleFile2, sep="\t", stringsAsFactors = F)
groups=split(groups_tbl$V2, groups_tbl$V1)
obj$group = unlist(groups[obj$orig.ident])

ngroup=length(unique(groups_tbl$V2))

genes_tbl<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
gene=genes_tbl$V1[1]
for (gene in genes_tbl$V1){
  gdata<-FetchData(object = obj, gene)
  colnames(gdata)<-"Gene"
  coords<-data.frame(obj@reductions$umap@cell.embeddings)
  gdata<-cbind(coords, data.frame(group=obj$group), gdata)
  gdata1<-subset(gdata, gdata$Gene == 0)
  gdata2<-subset(gdata, gdata$Gene > 0)
  
  png(filename=paste0(outFile, gene, ".feature_plot.png"), width= ngroup * 2000, height=2000, res=300)
  g<-ggplot(gdata, aes(UMAP_1, UMAP_2)) + 
    geom_point(data=gdata1, aes(col=Gene)) + 
    geom_point(data=gdata2, aes(col=Gene)) + 
    scale_colour_gradient(name=gene, low="grey", high="blue") + 
    facet_grid(~group) + 
    theme_bw() +
    theme(strip.background=element_rect(fill="white"))
  print(g)
  dev.off()
}
