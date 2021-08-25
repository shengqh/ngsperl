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
pca_dims<-1:npcs

prefix<-outFile

finalList<-readRDS(parFile1)
obj<-finalList$obj

obj<-run_cluster_only(obj, pca_dims, resolution, random.seed, reduction=reduction)

tbl<-table(obj$seurat_clusters, obj$orig.ident)
write.csv(tbl, file=paste0(outFile, ".cluster_sample.csv"))

seurat_clusters<-unlist(obj[["seurat_clusters"]])
seurat_colors<-hue_pal()(length(levels(seurat_clusters)))
names(seurat_colors)<-levels(seurat_clusters)

finalList$seurat_colors<-seurat_colors

counts=GetAssayData(obj,assay="RNA",slot="counts")
clusters=obj@active.ident
sumcounts=get_cluster_count(counts, clusters)
logsumcounts<-log2(sumcounts+1)

data.quantileAll<-apply(logsumcounts, 2, function(x){quantile(x, 0.75)})

norm_method=""
if(any(data.quantileAll == 0)){
  norm_method = ".normByTotal"
  data.all <- apply(logsumcounts, 2, sum)
  data.all<-data.all / median(data.all)
  data.norm <- t(t(logsumcounts) / data.all)
}else{
  norm_method = ".normByUpQuantile"
  data.quantileAll<-data.quantileAll / median(data.quantileAll)
  data.norm <- t(t(logsumcounts) / data.quantileAll)
}

write.csv(sumcounts, file=paste0(prefix, ".cluster.count.csv"))
write.csv(data.norm, file=paste0(prefix, ".cluster", norm_method, ".csv"))

clusters<-data.frame("cell" = c(1:length(obj$seurat_clusters)), "seurat_clusters"=as.numeric(as.character(obj$seurat_clusters)), stringsAsFactors = F)
rownames(clusters)<-names(obj$seurat_clusters)
write.csv(clusters, file=paste0(prefix, ".cluster.csv"))

#finalList<-readRDS(finalListFile)
#obj=finalList$obj

p<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_clusters") + guides(colour = guide_legend(override.aes = list(size = 3), ncol=1))
png(paste0(outFile, ".cluster.png"), width=3300, height=3000, res=300)
print(p)
dev.off()

finalList$obj<-obj
finalListFile<-paste0(prefix, ".final.rds")
saveRDS(finalList, file=finalListFile)
