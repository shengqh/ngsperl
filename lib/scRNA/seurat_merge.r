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
npcs<-as.numeric(myoptions$pca_dims)
pca_dims<-1:npcs

by_sctransform<-ifelse(myoptions$by_sctransform == "1", TRUE, FALSE)
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

obj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(obj, myoptions, prefix)
obj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

DefaultAssay(obj)<-"RNA"

if(by_sctransform){
  obj<-do_sctransform(obj, vars.to.regress=vars.to.regress, return.only.var.genes=FALSE)
  assay="SCT"
}else{
  assay="RNA"
}

#no matter if we will use sctransform, we need normalized RNA assay for visualization and cell type annotation
#data slot for featureplot, dotplot, cell type annotation and scale.data slot for heatmap
obj<-do_normalization(obj, selection.method="vst", nfeatures=3000, vars.to.regress=vars.to.regress, scale.all=FALSE, essential_genes=essential_genes)

DefaultAssay(obj)<-assay

cat("RunPCA ... \n")
obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)

# #https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# pcs <- find_number_of_reduction(obj, reduction="pca")
# final_pcs=max(pcs, as.numeric(myoptions$pca_dims))
# writeLines(paste0("pcs\t", final_pcs), con=paste0(outFile, ".pca.txt"))
# cat(paste0("recommend pcs=", final_pcs))

output_ElbowPlot(obj, outFile, "pca")

cat("run_umap ... \n")
obj <- RunUMAP(object = obj, dims=pca_dims, verbose = FALSE)

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

cat("output_integration_dimplot ... \n")
output_integration_dimplot(obj, outFile, FALSE, myoptions$qc_genes)

library('rmarkdown')
rmarkdown::render("seurat_data.rmd",output_file=paste0(outFile,".data.html"))

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
