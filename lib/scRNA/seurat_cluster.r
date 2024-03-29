rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge/result/AK6383.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_cluster_res0.5/result')

### Parameter setting end ###

source("scRNA_func.r")
library(dplyr)
library(Seurat)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(rmdformats)
library(DT)
library(data.table)
library(digest)
library(heatmap3)
library(cowplot)
library(scales)
library(stringr)
library(htmltools)
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

resolution=as.numeric(myoptions$resolution)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
reduction<-myoptions$reduction
npcs<-as.numeric(myoptions$pca_dims)

if(file.exists(parFile2)){
  ot<-read.table(parFile2, sep="\t")
  npcs<-as.numeric(ot$V2[1])
}

pca_dims<-1:npcs

prefix<-outFile

obj<-read_object(parFile1)

obj<-run_cluster_only(obj, pca_dims, resolution, random.seed, reduction=reduction)

tbl<-table(obj$seurat_clusters, obj$orig.ident)
write.csv(tbl, file=paste0(outFile, ".cluster_sample.csv"))

osn<-tbl / rowSums(tbl) * 5000
osm<-data.frame(osn)
colnames(osm)<-c("Cluster", "Sample", "Cell")
osm$Cluster<-paste0("Cluster ", osm$Cluster)
osm$Cluster<-factor(osm$Cluster, levels=unique(osm$Cluster))
osm_perc = osm %>% dplyr::group_by(Cluster) %>% dplyr::mutate(Percent = Cell/sum(Cell))

height=max(1600, min(5000, nrow(tbl) * 100 ))
width=max(1600, min(5000, ncol(tbl) * 100 ))
png(paste0(outFile, ".cluster_sample_percByCluster.png"), width=width, height=height, res=300)
g<-ggplot(osm_perc, aes(x=Sample,y=Cluster)) + geom_point(mapping = aes_string(size = "Percent")) + theme_bw() + xlab("") + ylab("") + theme(axis.text.x=element_text(angle = 90))
print(g)
dev.off()

ttbl<-t(tbl)
osn<-ttbl / rowSums(ttbl) * 100
osm<-data.frame(osn)
colnames(osm)<-c("Sample", "Cluster", "Percent")
osm$Cluster<-paste0("Cluster ", osm$Cluster)
osm$Cluster<-factor(osm$Cluster, levels=unique(osm$Cluster))

height=max(1600, min(5000, ncol(ttbl) * 100 ))
width=max(1600, min(5000, nrow(ttbl) * 100 ))
png(paste0(outFile, ".cluster_sample_percBySample.png"), width=width, height=height, res=300)
g<-ggplot(osm, aes(x=Sample,y=Cluster)) + geom_point(mapping = aes_string(size = "Percent")) + theme_bw() + xlab("") + ylab("") + theme(axis.text.x=element_text(angle = 90))
print(g)
dev.off()

seurat_clusters<-unlist(obj[["seurat_clusters"]])
seurat_colors<-hue_pal()(length(levels(seurat_clusters)))
names(seurat_colors)<-levels(seurat_clusters)

counts=GetAssayData(obj,assay="RNA",slot="counts")
clusters=obj$seurat_clusters
sumcounts=get_cluster_count(counts, clusters)

write.csv(sumcounts, file=paste0(prefix, ".cluster.count.csv"))

clusters<-data.frame("cell" = c(1:length(obj$seurat_clusters)), "seurat_clusters"=as.numeric(as.character(obj$seurat_clusters)), "orig.ident"=obj$orig.ident, stringsAsFactors = F)
rownames(clusters)<-names(obj$seurat_clusters)
write.csv(clusters, file=paste0(prefix, ".cluster.csv"))

p<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_clusters") + guides(colour = guide_legend(override.aes = list(size = 3), ncol=1))
png(paste0(outFile, ".cluster.png"), width=3300, height=3000, res=300)
print(p)
dev.off()

if(length(unique(obj$orig.ident)) > 1){
  p2<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="orig.ident")
  p3<-p+p2
  png(paste0(outFile, ".umap.sample_cell.png"), width=6600, height=3000, res=300)
  print(p3)
  dev.off()
}else{
  png(paste0(outFile, ".umap.sample_cell.png"), width=3300, height=3000, res=300)
  print(p)
  dev.off()
}

saveRDS(obj@meta.data, file=paste0(prefix, ".meta.rds"))
saveRDS(obj, file=paste0(prefix, ".final.rds"))

#we don't want to exponentiated the reads in "data" slot before taking average, so we have to copy it to "counts" slot first.
mean_exp<-get_seurat_average_expression(obj, cluster_name="seurat_clusters", assay="RNA")
write.csv(mean_exp, file=paste0(prefix, ".cluster.meanexp.csv"))
