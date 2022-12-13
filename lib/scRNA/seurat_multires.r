rm(list=ls()) 
outFile='PH_combine'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/scratch/cqs/shengq2/paula_hurley_projects/20221115_scRNA_7467_benign_hg38/seurat_sct_harmony/result/PH_combine.final.rds'
parFile2=''
parFile3=''


setwd('/scratch/cqs/shengq2/paula_hurley_projects/20221115_scRNA_7467_benign_hg38/seurat_sct_harmony_multires_01_call/result')

### Parameter setting end ###

source("scRNA_func.r")
library(dplyr)
library(Seurat)
library(knitr)
library(kableExtra)
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

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
by_integration<-ifelse(myoptions$by_integration == "0", FALSE, TRUE)
reduction<-myoptions$reduction
by_harmony<-reduction=="harmony"
assay=get_assay(by_sctransform, by_integration, by_harmony)

npcs<-as.numeric(myoptions$pca_dims)
species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype_str<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
curated_markers_file=myoptions$curated_markers_file

layer=ifelse(is.null(myoptions$layer), "Layer4", myoptions$layer)

bubblemap_file=myoptions$bubblemap_file

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = species,
                             curated_markers_file = curated_markers_file,
                             HLA_panglao5_file = HLA_panglao5_file,
                             layer=layer,
                             remove_subtype_str = remove_subtype_str,
                             combined_celltype_file = parSampleFile2)

tiers<-ctdef$tiers

cell_activity_database<-ctdef$cell_activity_database

celltype_map<-ctdef$celltype_map

combined_ct<-ctdef$combined_celltypes

prefix<-outFile

if(!exists("obj")){
  obj<-read_object(parFile1)
}

resolutions=c(0.5, 1.0, 1.5)

curreduction=ifelse(by_harmony, "harmony", "pca")

DefaultAssay(obj)<-assay

cat("FindNeighbors...\n")
obj <- FindNeighbors(object=obj, reduction=curreduction, dims=c(1:npcs), verbose=FALSE)

cat("FindClusters...\n")
obj <- FindClusters(object=obj, verbose=FALSE, random.seed=random.seed, resolution=resolutions)

res_prefix<-paste0(assay, "_snn_res")
multi_res<-colnames(obj@meta.data)[grepl(paste0("^", res_prefix, ".+\\d$"), colnames(obj@meta.data))]

if(identical(multi_res, character(0))){
  stop(paste0("cannot find resolution with prefix: ", res_prefix))
}

cat("Predict cell type...\n")
cur_res = multi_res[1]
for(cur_res in multi_res){
  cat("  ", cur_res, "\n")
  data.norm=get_seurat_average_expression(obj, cur_res)
  
  predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
  res_values=unlist(FetchData(obj, cur_res))
  
  #store the raw cell types
  new.cluster.ids<-names(predict_celltype$max_cta)
  names(new.cluster.ids) <- colnames(data.norm)
  
  cur_celltype<-paste0(cur_res, "_rawcelltype")
  obj[[cur_celltype]] = unlist(new.cluster.ids[res_values])
  
  #store mapped cell types (summerize sub celltypes, for example T cells)
  new.cluster.ids<-celltype_map[names(predict_celltype$max_cta)]
  names(new.cluster.ids) <- colnames(data.norm)
  
  cur_celltype<-paste0(cur_res, "_celltype")
  res_values=unlist(FetchData(obj, cur_res))
  obj[[cur_celltype]] = unlist(new.cluster.ids[res_values])
}

raw_cts<-paste0(multi_res, "_rawcelltype")
multi_cts<-paste0(multi_res, "_celltype")

res_df<-data.frame("resolution"=resolutions, "cluster"=multi_res, "celltype"=multi_cts, "rawcelltype"=raw_cts)
write.csv(res_df, file=paste0(outFile, ".resolutions.csv"), row.names=F)

#umaplist<-RunMultipleUMAP(obj, curreduction=curreduction, cur_pca_dims=c(1:npcs))
#obj<-RunUMAP(obj, reduction=curreduction, dims=c(1:npcs))

cat("Visualization...\n")
cind<-1
for(cind in c(1:nrow(res_df))){
  cur_res = res_df$cluster[cind]
  raw_celltype = res_df$rawcelltype[cind]
  cur_celltype = res_df$celltype[cind]
  cat(cur_celltype, "\n")
  
  meta<-obj@meta.data
  meta<-meta[order(meta[,cur_res]),]  
  seurat_celltype<-paste0(meta[,cur_res], ":", meta[,cur_celltype])
  seurat_celltype<-factor(seurat_celltype, levels=unique(seurat_celltype))
  sname = paste0("seurat_", cur_celltype)
  meta[,sname] = seurat_celltype

  seurat_rawcelltype<-paste0(meta[,cur_res], ":", meta[,raw_celltype])
  seurat_rawcelltype<-factor(seurat_rawcelltype, levels=unique(seurat_rawcelltype))
  s_rawname = paste0("seurat_", raw_celltype)
  meta[,s_rawname] = seurat_rawcelltype

  meta<-meta[colnames(obj),]
  obj@meta.data = meta

  umap_width=2300
  dot_width=4000

  # raw cell type
  g<-get_dim_plot_labelby(obj, label.by=raw_celltype, reduction="umap", legend.title="")
  png(paste0(prefix, ".", raw_celltype, ".umap.png"), width=umap_width, height=2000, res=300)
  print(g)
  dev.off()

  g<-get_bubble_plot(obj, cur_res=NA, raw_celltype, bubblemap_file, assay="RNA", orderby_cluster = TRUE)
  png(paste0(prefix, ".", raw_celltype, ".dot.png"), width=dot_width, height=get_dot_height(obj, raw_celltype), res=300)
  print(g)
  dev.off()

  g<-get_celltype_marker_bubble_plot( obj = obj, 
                                      group.by = raw_celltype, 
                                      cellType = cell_activity_database$cellType,
                                      weight = cell_activity_database$weight,
                                      n_markers = 5, 
                                      combined_ct=combined_ct)

  png(paste0(prefix, ".", raw_celltype, ".ct_markers.bubbleplot.png"), width=dot_width, height=get_dot_height(obj, raw_celltype), res=300)
  print(g)
  dev.off()

  # cell type, all clusters
  g<-get_dim_plot(obj, group.by=cur_res, label.by=sname, reduction="umap", legend.title="")
  png(paste0(prefix, ".", cur_celltype, ".seurat.umap.png"), width=umap_width + 1000, height=2000, res=300)
  print(g)
  dev.off()

  g<-get_bubble_plot(obj, cur_res=cur_res, cur_celltype, bubblemap_file, assay="RNA", orderby_cluster = FALSE)
  png(paste0(prefix, ".", cur_celltype, ".seurat.dot.png"), width=dot_width, height=get_dot_height(obj, cur_res), res=300)
  print(g)
  dev.off()

  for(pct in unique(unlist(obj[[cur_celltype]]))){
    cells=colnames(obj)[obj[[cur_celltype]] == pct]
    subobj=subset(obj, cells=cells)
    g<-get_dim_plot(subobj, group.by=cur_res, label.by=s_rawname, reduction="umap", legend.title="")

    png(paste0(prefix, ".", cur_celltype, ".", celltype_to_filename(pct), ".umap.png"), width=umap_width + 1000, height=2000, res=300)
    print(g)
    dev.off()
  }

  output_celltype_figures(obj, cur_celltype, prefix, bubblemap_file, cell_activity_database, combined_ct, group.by="orig.ident", name="sample")
  if("batch" %in% colnames(obj@meta.data)){
    output_celltype_figures(obj, cur_celltype, prefix, bubblemap_file, cell_activity_database, combined_ct, group.by="batch", name="batch")
  }

  if(!is.na(combined_ct)){
    new_layer = paste0(cur_celltype, "_summary")
    new_celltype<-as.character(obj@meta.data[,cur_celltype])
    in_new_celltype<-new_celltype %in% names(combined_ct)
    new_celltype[in_new_celltype]<-combined_ct[new_celltype[in_new_celltype]]
    obj[[new_layer]] = new_celltype

    output_celltype_figures(obj, new_layer, prefix, bubblemap_file, cell_activity_database, combined_ct, group.by="orig.ident", name="sample")
    if("batch" %in% colnames(obj@meta.data)){
      output_celltype_figures(obj, new_layer, prefix, bubblemap_file, cell_activity_database, combined_ct, group.by="batch", name="batch")
    }
  }
}

saveRDS(obj@meta.data, paste0(prefix, ".meta.rds"))
write.csv(obj@meta.data, paste0(prefix, ".meta.csv"))

library('rmarkdown')
rmarkdown::render("seurat_multires.rmd",output_file=paste0(outFile,".html"))

#saveRDS(obj, paste0(prefix, ".multires.rds"))
