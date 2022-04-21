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
library(glmGamPoi)
require(data.table)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "1", TRUE, FALSE)
sctransform_by_individual<-ifelse(myoptions$sctransform_by_individual == "1", TRUE, FALSE)
regress_by_percent_mt<-ifelse(myoptions$regress_by_percent_mt == "1", TRUE, FALSE)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile2, sep="\t" ,header=F)$V1

prefix<-outFile

species=myoptions$species

finalListFile<-paste0(prefix, ".final.rds")

rawobj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(rawobj, myoptions, prefix)
rawobj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

if(by_sctransform){
  cat("performing SCTransform ...\n")
  nsamples=length(unique(rawobj$sample))
  if(nsamples > 1 & sctransform_by_individual){
    objs<-SplitObject(object = rawobj, split.by = "sample")
    rm(rawobj)
  
    #perform sctransform
    objs<-lapply(objs, function(x){
      
      x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=TRUE, verbose = FALSE)
      return(x)
    })  
    obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
    VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data)
  }else{
    obj=rawobj
    rm(rawobj)
    obj<-SCTransform(obj, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=TRUE, verbose = FALSE)
  }
  assay="SCT"
}else{
  cat("performing NormalizeData/FindVariableFeatures ...\n")
  #perform standard workflow
  obj <-rawobj
  rm(rawobj)
  obj<-NormalizeData(obj, verbose = FALSE)
  obj<-FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE) 
  vgenes<-VariableFeatures(obj, selection.method = "vst")
  sgenes<-unique(c(essential_genes, vgenes))
  obj<-ScaleData(obj,vars.to.regress = vars.to.regress, features = sgenes)
  assay="RNA"
}

cat("run_pca ... \n")
obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)

#https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pcs <- find_number_of_reduction(obj, reduction="pca")
final_pcs=max(pcs, as.numeric(myoptions$pca_dims))
writeLines(paste0("pcs\t", final_pcs), con=paste0(outFile, ".pca.txt"))
cat(paste0("pcs=", final_pcs))

png(paste0(outFile, ".elbowplot.pca.png"), width=1500, height=1200, res=300)
p<-ElbowPlot(obj, ndims = 40, reduction = "pca")  + geom_vline(xintercept=final_pcs, color="red")
if (final_pcs != pcs){
  p<-p + geom_vline(xintercept=pcs, color="blue")
}
print(p)
dev.off()

pca_dims<-1:max(pcs, as.numeric(myoptions$pca_dims))

cat("run_umap ... \n")
obj <- RunUMAP(object = obj, dims=pca_dims, verbose = FALSE)

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

output_integration_dimplot(obj, outFile, FALSE)
