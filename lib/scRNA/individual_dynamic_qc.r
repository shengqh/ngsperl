rm(list=ls()) 
sample_name='Adipose_9240'
outFile='Adipose_9240'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1='/data/cqs/shengq2/program/collaborations/celestine_wanjalla/20230115_combined_scRNA_hg38/20230501_filter_config.txt'
parFile2=''
parFile3='/data/wanjalla_lab/projects/20231025_combined_scRNA_hg38_CITEseq/essential_genes/result/combined.txt'


setwd('/data/wanjalla_lab/projects/20231025_combined_scRNA_hg38_CITEseq/raw_dynamic_qc_sct2/result/Adipose_9240')

### Parameter setting end ###

source("scRNA_func.r")
library(plyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(digest)
library(cowplot)
library(scales)
library(stringr)
library(htmltools)
library(patchwork)
library(glmGamPoi)
library('rmarkdown')

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-is_one(myoptions$by_sctransform)
npcs<-as.numeric(myoptions$pca_dims)

species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype_str<-myoptions$remove_subtype
annotate_tcell<-is_one(myoptions$annotate_tcell)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(by_sctransform, "SCT", "RNA")
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)
resolution=as.numeric(myoptions$dynamic_by_one_resolution)
curated_markers_file=myoptions$curated_markers_file

Mtpattern=myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern

ensembl_map=NULL
if("ensembl_gene_map_file" %in% names(myoptions)){
  ensembl_gene_map_file = myoptions$ensembl_gene_map_file
  if(ensembl_gene_map_file != ""){
    gene_tb=read.table(ensembl_gene_map_file, sep="\t", header=T)
    gene_tb=gene_tb[!duplicated(gene_tb$ENSEMBL_GENE_ID),]
    gene_tb=gene_tb[gene_tb$ENSEMBL_GENE_ID != "",]
    ensembl_map = split(gene_tb$GENE_SYMBOL, gene_tb$ENSEMBL_GENE_ID)
  }
}

by_individual_sample=TRUE
by_harmony<-FALSE
redo_harmony<-FALSE

layer=ifelse(is.null(myoptions$layer), "Layer4", myoptions$layer)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is_file_empty(bubblemap_file)

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

if(file.exists(parFile3)){
  essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1
}else{
  essential_genes=NULL
}

if(!exists('parSampleFile4')){
  parSampleFile4=""
}

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = species,
                             curated_markers_file = curated_markers_file,
                             HLA_panglao5_file = HLA_panglao5_file,
                             layer=layer,
                             remove_subtype_str = remove_subtype_str,
                             combined_celltype_file = parSampleFile4)

cell_activity_database<-ctdef$cell_activity_database

layer2map<-ctdef$celltype_map

combined_ct<-ctdef$combined_celltypes
combined_ct_source<-ctdef$combined_celltype_source

replace_cts=layer2map %in% names(combined_ct)
layer2map[replace_cts] = combined_ct[layer2map[replace_cts]]

prefix<-outFile
root_folder=getwd()

raw_files=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
file_path=raw_files$V1[1]
sample_name=raw_files$V2[1]

obj = read_object_from_rawfile(sample_name, file_path, species, ensembl_map=NULL)

obj<-PercentageFeatureSet(object=obj, pattern=Mtpattern, col.name="percent.mt", assay="RNA")
obj<-PercentageFeatureSet(object=obj, pattern=rRNApattern, col.name = "percent.ribo", assay="RNA")

finalList= preprocessing_rawobj(obj, myoptions, prefix)
obj = finalList$rawobj

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes, species = myoptions$species)
  bubble_genes<-unique(genes_df$gene)
}

obj[["layer0"]]<-"Unassigned"
obj[["layer0_clusters"]]<-0
obj[["layer0_raw"]]<-"Unassigned"

cbind_celltype<-function(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts){
  if(is.null(cur_layermap)){
    return(cur_cts)
  }
  layer_ids<-unlist(cur_layermap[new_cluster_ids])
  names(layer_ids) <- colnames(data_norm)
  
  oldcluster<-subobj[[cluster]][[1]]
  cur_cts$seurat_clusters=oldcluster
  cur_cts$raw_cell_type<-new_cluster_ids[oldcluster]
  cur_cts$raw_seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$raw_cell_type) 
  cur_cts$cell_type<-layer_ids[oldcluster]
  cur_cts$seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$cell_type)

  return(cur_cts)
}

get_empty_files<-function(){
  files = data.frame("previous_layer"=character(),
                     "cur_layer"=character(),
                     "pct"=character(),
                     "type"=character(),
                     "fname"=character())
  return(files)
}

if(0){
  previous_celltypes<-unique(obj@meta.data[[previous_layer]])
  previous_layer<-"layer0"
  cur_layer="layer4"
  cur_layermap=layer2map
  iter=1
  iter_name=paste0("iter", iter)
  previous_celltypes<-previous_celltypes[order(previous_celltypes)]
  curprefix = paste0(prefix, ".iter", iter)
  tmp_folder = paste0(root_folder, "/details")
  if(!dir.exists(tmp_folder)){
    dir.create(tmp_folder)
  }
  setwd(tmp_folder)
}

if(0){
  previous_layer = "layer0"
  cur_layer = "layer4"
  cur_layermap = layer2map
}

if("harmony" %in% names(obj@reductions)){
  obj@reductions["harmony"]<-NULL
}
if("umap" %in% names(obj@reductions)){
  obj@reductions["umap"]<-NULL
}
if("pca" %in% names(obj@reductions)){
  obj@reductions["pca"]<-NULL
}

DefaultAssay(obj)<-"RNA"
if(by_sctransform){
  obj<-do_sctransform(obj, vars.to.regress=vars.to.regress, return.only.var.genes=FALSE, mc.cores=8, use_sctransform_v2=TRUE)
}

assay=ifelse(by_sctransform, "SCT", "RNA")

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

cat("RunUMAP ... \n")
obj <- RunUMAP(object = obj, dims=pca_dims, verbose = FALSE)

if(by_sctransform){
  #clear SCT data for small object
  #https://github.com/satijalab/seurat/issues/2587
  obj[["SCT"]]@misc <- NULL    
}    

cur_folder = getwd()
tmp_folder = paste0(cur_folder, "/details")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}
obj <- do_analysis( tmp_folder = tmp_folder,
                    cur_folder = cur_folder,
                    obj = obj, 
                    layer2map = layer2map, 
                    npcs = npcs, 
                    resolution = resolution, 
                    random.seed = random.seed, 
                    by_sctransform = by_sctransform, 
                    by_harmony = by_harmony, 
                    prefix = prefix, 
                    vars.to.regress = vars.to.regress, 
                    bubblemap_file = bubblemap_file, 
                    essential_genes = essential_genes,
                    by_individual_sample = FALSE,
                    species = species,
                    reduction="pca")

saveRDS(obj, file=paste0(prefix, ".obj.rds"))

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
