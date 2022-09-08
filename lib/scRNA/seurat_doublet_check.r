
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

if(!exists('obj')){
  obj=readRDS(parFile1)
  if(is.list(obj)){
    obj<-obj$obj
  }
  meta1=readRDS(parFile2)
  meta1$id=rownames(meta1)
  meta2=readRDS(parFile3)
  meta2$id=rownames(meta2)
  meta=merge(meta1, meta2, by="id")
  rownames(meta)<-meta$id
  meta=meta[colnames(obj),]
  obj@meta.data=meta
}

doublet_options=read.csv(parFile4)

all_dt=NULL
idx=1
for(idx in c(1:nrow(doublet_options))){
  doublet_rate=doublet_options$doublet_rate[idx]
  doublet_col=doublet_options$label[idx]

  dt=as.data.frame.matrix(table(obj@meta.data[,myoptions$cluster_celltype_layer], obj@meta.data[,doublet_col]))
  dt$Total<-dt$Doublet+dt$Singlet
  dt$DoubletPerc=dt$Doublet / dt$Total
  dt$DoubletRate=doublet_rate
  dt$Cluster=rownames(dt)

  all_dt<-rbind(all_dt, dt)
  
  doublet_cells=colnames(obj)[obj@meta.data[,doublet_col]=="Doublet"]
  
  g1<-DimPlot(obj, group.by = myoptions$cluster_layer, label=T) + ggtitle(myoptions$celltype_layer)+
    scale_color_discrete(labels = levels(obj@meta.data[,myoptions$cluster_celltype_layer]))
  
  g2<-DimPlot(obj, label=F, cells.highlight =doublet_cells) + 
    ggtitle("DoubletFinder") + scale_color_discrete(type=c("gray", "red"), labels = c("Singlet", "Doublet"))
  
  if(has_bubblemap){
    subobj=subset(obj, cells=doublet_cells)
    layout <- "
AABB
CCCC
"
    g3<-get_bubble_plot(subobj, myoptions$cluster_layer, myoptions$celltype_layer, bubblemap_file, assay="RNA")
    g<-g1+g2+g3+plot_layout(design=layout)
    height=4000
  }else{
    g<-g1+g2+plot_layout(ncol=2)
    height=2000
  }
  
  png(paste0(outFile, ".DR", sprintf(doublet_rate, fmt="%#.2f"), ".umap.png"), width=5400, height=height, res=300)
  print(g)
  dev.off()
}

all_dt$Cluster<-factor(all_dt$Cluster, levels=levels(obj@meta.data[,myoptions$cluster_celltype_layer]))
png(paste0(outFile, ".doublet_perc.png"), width=3000, height=3000, res=300)
g<-ggplot(all_dt, aes(x=DoubletRate, y=DoubletPerc)) + geom_line()  + facet_wrap(~Cluster) + theme_bw3()
print(g)
dev.off()

write.csv(all_dt, paste0(outFile, ".doublet_perc.csv"))

