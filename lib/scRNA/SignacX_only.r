rm(list=ls()) 
sample_name='Aorta_9240'
outFile='Aorta_9240'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/raw_qc_SignacX/result/Aorta_9240')

### Parameter setting end ###

source("scRNA_func.r")
library(SignacX)
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

options(Seurat.object.assay.version = 'v3')
options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
pca_dims=1:as.numeric(myoptions$pca_dims)
reduction=myoptions$reduction
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

find_neighbors <- is_one(myoptions$find_neighbors)

assay=myoptions$assay
if(assay == "RNA"){
  assay=ifelse(by_sctransform, "SCT", "RNA")  
}

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

DefaultAssay(obj) <- assay

if(is_seurat_5_plus(obj)){
  if(DefaultAssay(obj) == "integrated"){
    if(nrow(obj@assays$integrated@counts) == 0){
      if ("SCT" %in% names(obj@assays)){
        counts = MyGetAssayData(obj, "SCT", slot="counts")
      }else{
        counts = MyGetAssayData(obj, "RNA", slot="counts")
      }
    }else{
      counts = MyGetAssayData(obj, "integrated", slot="counts")
    }
  }else{
    counts = MyGetAssayData(obj, assay, slot="counts")
  }
  newobj=CreateSeuratObject(counts, assay="RNA")
  newobj@reductions <- obj@reductions
  newobj<-FindNeighbors(object = newobj, reduction=reduction, dims=pca_dims, verbose=FALSE)
  newobj@meta.data=obj@meta.data

  obj = newobj
  assay = "RNA"
  rm(newobj)
}

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

ggsave(paste0(outFile, ".SignacX.png"), g, width=width, height=height, units="px", dpi=300, bg="white")

ct<-data.frame("SignacX"=obj$signacx_CellStates, "Sample"=obj$orig.ident)
ct_tbl<-table(ct$SignacX,ct$Sample)
write.csv(ct_tbl, paste0(outFile, ".SignacX_Sample.csv"))

unlink('.cache', recursive = TRUE, force = TRUE)

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))
