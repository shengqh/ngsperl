rm(list=ls()) 
sample_name='cvd_10a'
outFile='cvd_10a'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250211_T04_snRNA_hg38/raw_qc_sct2_SignacX/result/cvd_10a')

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

bubblemap_width=to_numeric(myoptions$bubblemap_width, 4000)
bubblemap_height=to_numeric(myoptions$bubblemap_height, 2000)
bubblemap_unit=ifelse(bubblemap_width > 50, "px", "in")

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
  obj=read_object_from_file_list(file_list_path=parSampleFile1)
}

DefaultAssay(obj) <- assay

# if(is_seurat_5_plus(obj)){
#   if(DefaultAssay(obj) == "integrated"){
#     if(nrow(obj@assays$integrated@counts) == 0){
#       if ("SCT" %in% names(obj@assays)){
#         counts = MyGetAssayData(obj, "SCT", slot="counts")
#       }else{
#         counts = MyGetAssayData(obj, "RNA", slot="counts")
#       }
#     }else{
#       counts = MyGetAssayData(obj, "integrated", slot="counts")
#     }
#   }else{
#     counts = MyGetAssayData(obj, assay, slot="counts")
#   }
#   newobj=CreateSeuratObject(counts, assay="RNA")
#   newobj@reductions <- obj@reductions
#   newobj<-FindNeighbors(object = newobj, reduction=reduction, dims=pca_dims, verbose=FALSE)
#   newobj@meta.data=obj@meta.data

#   obj = newobj
#   assay = "RNA"
#   rm(newobj)
# }

labels <- Signac(E=obj, R=R)

celltypes = GenerateLabels(labels, E = obj)
saveRDS(celltypes, file=paste0(outFile, ".SignacX.rds"))

ct_name="signacx_CellStates"
obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = ct_name)

ct<-data.frame("SignacX"=obj$signacx_CellStates, "Sample"=obj$orig.ident)
ct_tbl<-table(ct$SignacX,ct$Sample)
write.csv(ct_tbl, paste0(outFile, ".SignacX_Sample.csv"))

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

major_obj=get_category_with_min_percentage(obj, ct_name, 0.01)
ct_name_count = paste0(ct_name, "_count")
major_obj@meta.data = add_column_count(major_obj@meta.data, ct_name, ct_name_count)

g=get_dim_plot_labelby(major_obj, label.by = ct_name, reduction="umap", pt.size=0.1) + theme(plot.title=element_blank())
ggsave(paste0(outFile, ".SignacX.png"), g, width=6, height=4, units="in", dpi=300, bg="white")

if(has_bubblemap){
  g<-get_bubble_plot(
    obj=major_obj, 
    cur_res=NA, 
    cur_celltype=ct_name_count, 
    bubblemap_file, 
    assay="RNA", 
    species=myoptions$species,
    dot.scale=4)
    
  ggsave( paste0(outFile, ".SignacX.dot.png"), 
          g, 
          width=bubblemap_width, 
          height=bubblemap_height, 
          units=bubblemap_unit, 
          dpi=300, 
          bg="white")
}
rm(major_obj)

unlink('.cache', recursive = TRUE, force = TRUE)

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))
