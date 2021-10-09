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

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
prefix<-outFile

species=myoptions$species

pca_dims<-1:as.numeric(myoptions$pca_dims)

finalListFile<-paste0(prefix, ".final.rds")

rawobj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(rawobj, myoptions, prefix)
rawobj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

if(by_sctransform){
  cat("performing SCTransform ...\n")
  nsamples=length(unique(rawobj$sample))
  if(nsamples > 1){
    objs<-SplitObject(object = rawobj, split.by = "sample")
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
  obj<-ScaleData(obj)
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

output_integraion_dimplot(obj, outFile, FALSE)
