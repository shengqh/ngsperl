rm(list=ls()) 
outFile='pbmc_rejection'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1='/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250226_T04_snRNA_hg38/seurat_sct2_fastmnn/result/pbmc_rejection.final.rds'
parFile2=''
parFile3='/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250226_T04_snRNA_hg38/essential_genes/result/pbmc_rejection.txt'


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250226_T04_snRNA_hg38/seurat_sct2_fastmnn_dr0.5_1_call/result')

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
by_harmony<-reduction=="harmony"
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)
redo_harmony<-is_one(myoptions$redo_harmony, 0)
resolution=as.numeric(myoptions$dynamic_by_one_resolution)
curated_markers_file=myoptions$curated_markers_file
by_individual_sample=is_one(myoptions$by_individual_sample)

reduction=myoptions$reduction

if(by_individual_sample){
  by_harmony<-FALSE
}

layer=ifelse(is.null(myoptions$layer), "Layer4", myoptions$layer)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is_file_empty(bubblemap_file)

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
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
  obj<-readRDS(parFile1)
  if(is.list(obj)){
    obj<-obj$obj
  }
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


if(by_individual_sample){
  if("harmony" %in% names(obj@reductions)){
    obj@reductions["harmony"]<-NULL
  }
  if("umap" %in% names(obj@reductions)){
    obj@reductions["umap"]<-NULL
  }
  if("pca" %in% names(obj@reductions)){
    obj@reductions["pca"]<-NULL
  }
  samples <- unique(obj$orig.ident)
  result_list = c()
  all_ct_counts = NULL

  sample=samples[2]
  for (sample in samples){
    cat("processing ", sample, "...\n")
    cur_folder = paste0(root_folder, "/", sample)
    if(!dir.exists(cur_folder)){
      dir.create(cur_folder)
    }
    tmp_folder = paste0(cur_folder, "/details")
    if(!dir.exists(tmp_folder)){
      dir.create(tmp_folder)
    }
    
    f1<-read.table(paste0(root_folder, "/fileList1.txt"), sep="\t")
    f1$V1[f1$V2=="task_name"] = sample
    write.table(f1, paste0(cur_folder, "/fileList1.txt"), sep="\t", row.names=F)

    for(idx in c(2:9)){
      f=paste0(root_folder, "/fileList", idx, ".txt")
      if(file.exists(f)){
        file.copy(f, paste0(cur_folder, "/", basename(f)), overwrite = TRUE)
      }
    }

    file.copy(paste0(root_folder, "/scRNA_func.r"), paste0(cur_folder, "/scRNA_func.r"), overwrite = TRUE)
    file.copy(paste0(root_folder, "/seurat_scDynamic_one_layer_one_resolution.rmd"), paste0(cur_folder, "/seurat_scDynamic_one_layer_one_resolution.rmd"), overwrite = TRUE)
    file.copy(paste0(root_folder, "/reportFunctions.R"), paste0(cur_folder, "/reportFunctions.R"), overwrite = TRUE)

    subobj<-subset(obj, cells=colnames(obj)[obj@meta.data$orig.ident==sample])

    #remove genes without count
    sub_counts<-subobj@assays$RNA@counts
    rc<-rowSums(sub_counts)
    rc<-rc[rc>0]
    rm(sub_counts)
    subobj <- subset(subobj, features = names(rc))

    if(0){#for test
      subobj=subset(subobj, cells=colnames(subobj)[1:2000])
    }

    setwd(cur_folder)

    #for each sample, do its own PCA, FindClusters and UMAP first
    subobj = sub_cluster( subobj = subobj, 
                          assay = assay, 
                          by_sctransform = by_sctransform, 
                          by_harmony = FALSE, 
                          redo_harmony = FALSE, 
                          curreduction = "pca", 
                          k_n_neighbors = min(npcs, 20),
                          u_n_neighbors = min(npcs, 30),
                          random.seed = random.seed,
                          resolutions = resolution,
                          cur_npcs = npcs, 
                          cur_pca_dims = c(1:npcs),
                          vars.to.regress = vars.to.regress,
                          essential_genes = NULL)    

    res_list = do_analysis( tmp_folder = tmp_folder,
                            cur_folder = cur_folder,
                            obj = subobj, 
                            layer2map = layer2map, 
                            npcs = npcs, 
                            resolution = resolution, 
                            random.seed = random.seed, 
                            by_sctransform = by_sctransform, 
                            by_harmony = 0, 
                            prefix = sample, 
                            vars.to.regress = vars.to.regress, 
                            bubblemap_file = bubblemap_file, 
                            essential_genes = NULL,
                            by_individual_sample = 1,
                            species = species,
                            reduction=reduction)

    result_list<-c(result_list, res_list$html)
    all_ct_counts<-rbind(all_ct_counts, res_list$ct_count)
  }
  setwd(root_folder)
  writeLines(result_list, paste0(prefix, ".samples.list"))
  write.table(all_ct_counts, paste0(prefix, ".ct_count.tsv"), sep="\t", row.names=T)

  # output_file=paste0(prefix,".dynamic_individual.html")
  # rmdfile = "seurat_scDynamic_one_layer_one_resolution_summary.rmd"
  # rmarkdown::render(rmdfile, output_file=output_file)
}else{
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
                          essential_genes = NULL,
                          by_individual_sample = 0,
                          species = species,
                          reduction=reduction);
}
