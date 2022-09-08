source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(kableExtra)
library(dplyr)

options_table<-read.table("fileList1.txt", sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

finalList=readRDS(parFile1)
all_obj=finalList$obj
seurat_colors=finalList$seurat_colors

celltype=read.csv(parFile2)
celltype$seurat_cellactivity_clusters=paste0(celltype$seurat_clusters, " : ", celltype$cellactivity_clusters)
ctmap=split(celltype$seurat_cellactivity_clusters, celltype$seurat_clusters)

all_obj$seurat_cellactivity_clusters=unlist(ctmap[as.character(all_obj$seurat_clusters)])

gmt<-read.table(parFile3, sep=",", fill=NA, header=F)

idx=1
for(idx in 1:nrow(gmt)){
  line=gmt[idx,1]
  parts = unlist(strsplit(line, "\t"))
  name=parts[1]
  genes=parts[3:length(parts)]
  genes<-genes[genes %in% rownames(all_obj)]
  width=ceiling(length(genes) * 0.3)
  dot_filename=paste0(name, ".dot.pdf")
  pdf(file=dot_filename, width=width, height=7)
  g=DotPlot(all_obj, features=genes, assay="RNA", group.by="seurat_cellactivity_clusters" ) + 
    xlab("") + ylab("") + ggtitle(name) + 
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust=1))
  print(g)
  dev.off()
}
