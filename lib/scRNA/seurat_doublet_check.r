source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(scales)
library(stringr)
library(htmltools)
library(patchwork)

options(future.globals.maxSize= 10779361280)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

by_sctransform<-ifelse(myoptions$by_sctransform == "1", TRUE, FALSE)
assay=ifelse(by_sctransform, "SCT", "RNA")

prefix<-outFile

obj=readRDS(parFile1)
meta=readRDS(parFile2)
obj@meta.data=meta

doublet=read.csv(parFile3, row.names=1)
colname=colnames(doublet)[grepl('DF.classification', colnames(doublet))][1]
doublet=doublet[colnames(obj),]
obj=AddMetaData(obj, doublet[,colname], col.name="DF.classification")

dt=table(obj@meta.data[,myoptions$cluster_celltype_layer], obj$DF.classification)
write.csv(dt, paste0(outFile, ".doublet.csv"))

doublet_cells=colnames(obj)[obj$DF.classification=="Doublet"]

g1<-DimPlot(obj, group.by = myoptions$cluster_layer, label=T) + ggtitle(myoptions$celltype_layer)+
  scale_color_discrete(labels = levels(obj@meta.data[,myoptions$celltype_layer]))
g2<-DimPlot(obj, label=F, cells.highlight =doublet_cells) + 
  ggtitle("DoubletFinder") + scale_color_discrete(type=c("gray", "red"), labels = c("Singlet", "Doublet"))
g<-g1+g2+plot_layout(ncol=2)

png(paste0(outFile, ".doublet.umap.png"), width=4400, height=2000, res=300)
print(g)
dev.off()


if(has_bubblemap){
  g1<-get_bubble_plot(obj, myoptions$cluster_layer, myoptions$celltype_layer, bubblemap_file, assay)

  singlet_obj=subset(obj, DF.classification=="Singlet")
  g2<-get_bubble_plot(singlet_obj, myoptions$cluster_layer, myoptions$celltype_layer, bubblemap_file, assay)

  g<-g1+g2+plot_layout(ncol=1)
  png(paste0(outFile, ".bubble.png"), width=4400, height=4000, res=300)
  print(g)
  dev.off()
}
