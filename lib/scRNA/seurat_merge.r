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
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_rRNA<-ifelse(myoptions$Remove_rRNA == "0", FALSE, TRUE)
Remove_MtRNA<-ifelse(myoptions$Remove_MtRNA == "0", FALSE, TRUE)
resolution=as.numeric(myoptions$resolution)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
prefix<-outFile

nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
mt_cutoff=as.numeric(myoptions$mt_cutoff)
species=myoptions$species

pca_dims<-1:as.numeric(myoptions$pca_dims)

finalList<-list()
finalListFile<-paste0(prefix, ".final.rds")

rawobj<-readRDS(parFile1)

rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
Mt.genes<- grep (pattern= Mtpattern,rownames(rawobj), value=TRUE )

#filter cells
finalList$filter<-list(nFeature_cutoff_min=nFeature_cutoff_min,
                       nFeature_cutoff_max=nFeature_cutoff_max,
                       mt_cutoff=mt_cutoff,
                       nCount_cutoff=nCount_cutoff)

rawobj<-subset(rawobj, subset = nFeature_RNA > nFeature_cutoff_min & nFeature_RNA<nFeature_cutoff_max & nCount_RNA > nCount_cutoff & percent.mt < mt_cutoff)

if(Remove_rRNA){
  rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]
}

if(Remove_MtRNA){
  rawobj<-rawobj[!(rownames(rawobj) %in% Mt.genes),]
}

nsamples=length(unique(rawobj$orig.ident))

if(by_sctransform){
  cat("performing SCTransform ...\n")
  if(nsamples > 1){
    objs<-SplitObject(object = rawobj, split.by = "orig.ident")
    rm(rawobj)
  
    #perform sctransform
    objs<-lapply(objs, function(x){
      x <- SCTransform(x, verbose = FALSE)
      return(x)
    })  
    obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
    VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data)
  }else{
    obj=rawobj
    rm(rawobj)
    obj<-SCTransform(obj, verbose = FALSE)
  }
  assay="SCT"
}else{
  cat("performing NormalizeData/FindVariableFeatures ...\n")
  #perform standard workflow
  obj <-rawobj
  rm(rawobj)
  obj<-NormalizeData(obj, verbose = FALSE)
  obj<-FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)  
  assay="RNA"
}

cat("run_pca ... \n")
obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)

png(paste0(outFile, ".elbowplot.pca.png"), width=1500, height=1200, res=300)
p<-ElbowPlot(obj, ndims = 20, reduction = "pca")
print(p)
dev.off()

cat("run_umap ... \n")
obj <- RunUMAP(object = obj, dims=pca_dims, verbose = FALSE)

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

mt<-data.frame(UMAP_1=obj@reductions$umap@cell.embeddings[,1], 
               UMAP_2=obj@reductions$umap@cell.embeddings[,2],
               Sample=obj$orig.ident)

nSamples = length(unique(obj$orig.ident))
nWidth=ceiling(sqrt(nSamples))
nHeight=ceiling(nSamples / nWidth)

cat("draw pictures ... \n")
draw_dimplot(mt, paste0(outFile, ".merge_sample.png"), "Sample")

png(paste0(outFile, ".merge.png"), width=1500, height=1200, res=300)
p<-DimPlot(object = obj, reduction = 'umap', label=FALSE, group.by="orig.ident")
print(p)
dev.off()
