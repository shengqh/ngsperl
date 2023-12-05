rm(list=ls()) 
sample_name='DB_1'
outFile='DB_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20231202_scRNA_5798_human_liver_redo_cellbender/raw_qc_sct2_SignacX/result/DB_1')

### Parameter setting end ###

source("scRNA_func.r")
library(SignacX)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
pca_dims=1:as.numeric(myoptions$pca_dims)
reduction=myoptions$reduction
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
assay=ifelse(by_sctransform, "SCT", "RNA")

SignacX_reference_file=myoptions$SignacX_reference_file
if(is.null(SignacX_reference_file)){
  R="default"
}else if(SignacX_reference_file == ""){
  R="default"
}else if(!file.exists(SignacX_reference_file)){
  stop(paste0('reference file not exists: ', SignacX_reference_file))
}else{
  R = readRDS(SignacX_reference_file)
}

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

if(DefaultAssay(obj) == "integrated"){
  if(nrow(obj@assays$integrated@counts) == 0){
    if ("SCT" %in% names(obj@assays)){
      obj@assays$integrated@counts = obj@assays$SCT@counts
    }else{
      obj@assays$integrated@counts = obj@assays$RNA@counts
    }
  }
}

obj<-FindNeighbors(object = obj, reduction=reduction, dims=pca_dims, verbose=FALSE)

labels <- Signac(E=obj, R=R)

celltypes = GenerateLabels(labels, E = obj)
saveRDS(celltypes, file=paste0(outFile, ".SignacX.rds"))

obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "signacx_CellStates")

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

g1=MyDimPlot(obj, group.by = "signacx_CellStates", reduction="umap", label=T)

if(has_bubblemap){
  g2<-get_bubble_plot(obj, NA, "signacx_CellStates", bubblemap_file, assay="RNA", 
    species=myoptions$species)
  layout <- "ABB"
  g<-g1+g2+plot_layout(design=layout)
  width=6300
}else{
  g<-g1
  width=2300
}
height=2000

png(paste0(outFile, ".SignacX.png"), width=width, height=height, res=300)
print(g)
dev.off()

ct<-data.frame("SignacX"=obj$signacx_CellStates, "Sample"=obj$orig.ident)
ct_tbl<-table(ct$SignacX,ct$Sample)
write.csv(ct_tbl, paste0(outFile, ".SignacX_Sample.csv"))

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))
