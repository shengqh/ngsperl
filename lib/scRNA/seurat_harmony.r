source("scRNA_func.r")

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
npcs<-as.numeric(myoptions$pca_dims)
pca_dims<-1:npcs

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

obj=do_harmony(objs, by_sctransform, Remove_Mt_rRNA, Mtpattern, rRNApattern, npcs, parSampleFile2)
reduction="harmony"
rm(objs)

for (reduct in c("pca", "harmony")){
  png(paste0(outFile, ".elbowplot.", reduct, ".png"), width=1500, height=1200, res=300)
  p<-ElbowPlot(obj, ndims = 20, reduction = reduct)
  print(p)
  dev.off()
}

cat("run_umap ... ")
obj <- RunUMAP(object = obj, reduction=reduction, dims=pca_dims, verbose = FALSE)

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

mt<-data.frame(UMAP_1=obj@reductions$umap@cell.embeddings[,1], 
               UMAP_2=obj@reductions$umap@cell.embeddings[,2],
               Sample=obj$orig.ident,
               batch=obj$batch)

nSamples = length(unique(obj$orig.ident))
nWidth=ceiling(sqrt(nSamples))
nHeight=ceiling(nSamples / nWidth)

cat("draw pictures ... ")
draw_dimplot(mt, paste0(outFile, ".harmony_sample.png"), "Sample")
if(has_batch_file){
  draw_dimplot(mt, paste0(outFile, ".harmony_batch.png"), "batch")
}

p1<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="orig.ident")
width=1500
if(has_batch_file){
  p2<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="batch")
  p<-p1+p2
  width=3000
}
png(paste0(outFile, ".harmony.png"), width=width, height=1200, res=300)
print(p)
dev.off()
