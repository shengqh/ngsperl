source('scRNA_func.r')

library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(cowplot)
library(scales)
library(stringr)
library(harmony)
library(patchwork)
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

optionstr<-paste0(options_table$V1, collapse = "_")

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_rRNA<-ifelse(myoptions$Remove_rRNA == "0", FALSE, TRUE)
Remove_MtRNA<-ifelse(myoptions$Remove_MtRNA == "0", FALSE, TRUE)
resolution=as.numeric(myoptions$resolution)
batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

prefix<-outFile

nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
mt_cutoff=as.numeric(myoptions$mt_cutoff)
species=myoptions$species
#nCount_sd_cutoff=as.numeric(options$nCount_sd_cutoff)
nCount_sd_cutoff=0

has_batch_file<-file.exists(parSampleFile2)

pca_dims<-1:as.numeric(myoptions$pca_dims)

finalList<-list()
finalListFile<-paste0(prefix, ".final.rds")

rawobj<-readRDS(parFile1)

rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
Mt.genes<- grep (pattern= Mtpattern,rownames(rawobj), value=TRUE )

rawobj<-subset(rawobj, subset = nFeature_RNA > nFeature_cutoff_min & nFeature_RNA<nFeature_cutoff_max & nCount_RNA > nCount_cutoff & percent.mt < mt_cutoff)

if(Remove_rRNA){
  rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]
}

if(Remove_MtRNA){
  rawobj<-rawobj[!(rownames(rawobj) %in% Mt.genes),]
}

objs<-SplitObject(object = rawobj, split.by = "orig.ident")
rm(rawobj)

#filter cells
finalList$filter<-list(nFeature_cutoff_min=nFeature_cutoff_min,
                      nFeature_cutoff_max=nFeature_cutoff_max,
                      mt_cutoff=mt_cutoff,
                      nCount_cutoff=nCount_cutoff,
                      nCount_sd_cutoff=nCount_sd_cutoff)

obj=do_harmony(objs, by_sctransform, parSampleFile2)
reduction="harmony"
rm(objs)

cat("run_cluster_only ... ")
obj <- run_cluster_only(obj, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, reduction=reduction)

cat("saving result ... ")
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

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

mt<-data.frame(UMAP_1=obj@reductions$umap@cell.embeddings[,1], 
               UMAP_2=obj@reductions$umap@cell.embeddings[,2],
               Sample=obj$orig.ident,
               batch=obj$batch)

if(has_batch_file){
  mt$batch=obj$batch
}

nSamples = length(unique(obj$orig.ident))
nWidth=ceiling(sqrt(nSamples))
nHeight=ceiling(nSamples / nWidth)

cat("draw pictures ... ")
draw_dimplot(mt, paste0(outFile, ".cluster_sample.png"), "Sample")
draw_dimplot(mt, paste0(outFile, ".cluster_batch.png"), "batch")

p1<-DimPlot(object = obj, reduction = 'umap', label=TRUE, group.by="seurat_clusters", cols=seurat_colors)
p2<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="orig.ident")
p<-p1+p2
width=3000
if(has_batch_file){
  p3<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="batch")
  p<-p1+p3+p2
  width=4500
}
png(paste0(outFile, ".cluster.png"), width=width, height=1200, res=300)
print(p)
dev.off()
