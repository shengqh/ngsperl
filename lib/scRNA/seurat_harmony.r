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

prefix<-outFile

has_batch_file<-file.exists(parSampleFile2)
npcs<-as.numeric(myoptions$pca_dims)
pca_dims<-1:npcs

obj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(obj, myoptions, prefix)
obj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

essential_genes=read.table(parFile2, sep="\t" ,header=F)$V1

by_sctransform<-ifelse(myoptions$by_sctransform == "1", TRUE, FALSE)
regress_by_percent_mt<-ifelse(myoptions$regress_by_percent_mt == "1", TRUE, FALSE)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

obj<-do_harmony(obj, by_sctransform, vars.to.regress, has_batch_file, parSampleFile2, pca_dims, essential_genes=essential_genes, mc.cores=8)

reduction="harmony"

for (reduct in c("pca", "harmony")){
  output_ElbowPlot(obj, outFile, reduct)
}

cat("RunUMAP ... ")
obj <- RunUMAP(object = obj, reduction=reduction, dims=pca_dims, verbose = FALSE)

finalList$obj<-obj

finalListFile<-paste0(outFile, ".final.rds")
saveRDS(finalList, file=finalListFile)

output_integration_dimplot(obj, outFile, FALSE, myoptions$qc_genes)

