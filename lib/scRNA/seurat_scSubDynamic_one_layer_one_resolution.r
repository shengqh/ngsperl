rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1='/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/seurat_sct2_merge/result/iSGS_cell_atlas.final.rds'
parFile2='/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.1_1_call/result/iSGS_cell_atlas.scDynamic.meta.rds'
parFile3='/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/essential_genes/result/iSGS_cell_atlas.txt'


setwd('/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.1_2_call_sub/result')

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
library(heatmap3)
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

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-is_one(myoptions$by_sctransform)
reduction<-myoptions$reduction
npcs<-as.numeric(myoptions$pca_dims)

species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype_str<-myoptions$remove_subtype
annotate_tcell<-is_one(myoptions$annotate_tcell)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(by_sctransform, "SCT", "RNA")
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)
by_harmony<-is_one(myoptions$redo_harmony)
resolution=as.numeric(myoptions$dynamic_by_one_resolution)
curated_markers_file=myoptions$curated_markers_file

previous_layer=myoptions$init_layer
cur_layer=myoptions$final_layer

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is_file_empty(bubblemap_file)

pca_dims<-1:npcs

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

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2)
  Idents(obj)<-init_layer
}

if(parSampleFile2 != ""){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  cat("removing genes in", ignore_gene_files$V1, "\n")
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  obj<-obj[!(rownames(obj) %in% ignore_genes),]
}

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes, species = myoptions$species)
  bubble_genes<-unique(genes_df$gene)
}

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
  cur_layermap=layer2map
  iter=1
  iter_name=paste0("sub_iter", iter)
  previous_celltypes<-previous_celltypes[order(previous_celltypes)]
  curprefix = paste0(prefix, ".sub_iter", iter)
  tmp_folder = paste0(root_folder, "/details")
  if(!dir.exists(tmp_folder)){
    dir.create(tmp_folder)
  }
  setwd(tmp_folder)
}

cur_folder = getwd()
tmp_folder = paste0(cur_folder, "/details")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}
res_list <- do_analysis(tmp_folder = tmp_folder,
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
                        by_individual_sample = 0,
                        species = species,
                        init_layer=previous_layer,
                        final_layer=cur_layer);

