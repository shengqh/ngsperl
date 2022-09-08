rm(list=ls()) 
outFile='AG_integrated'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='C:/projects/data/h_gelbard_lab/projects/20220803_integrated_project/seurat_sct_harmony/result/AG_integrated.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/data/h_gelbard_lab/projects/20220803_integrated_project/seurat_sct_harmony_multires_01_call/result')

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
reduction<-myoptions$reduction
npcs<-as.numeric(myoptions$pca_dims)
species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype_str<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
curated_markers_file=myoptions$curated_markers_file
assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")
by_harmony<-reduction=="harmony"

bubblemap_file=myoptions$bubblemap_file

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

if(file.exists(parSampleFile2)){
  cts<-read.table(parSampleFile2, header=F, sep="\t", stringsAsFactors = F)
  combined_ct<-unlist(split(cts$V2, cts$V1))
}

tiers<-read.table(myoptions$HLA_panglao5_file, sep="\t", header=T)

cell_activity_database<-read_cell_markers_file(
  panglao5_file = markerfile, 
  species = species, 
  remove_subtype_str = remove_subtype_str, 
  HLA_panglao5_file = HLA_panglao5_file, 
  curated_markers_file=curated_markers_file, 
  remove_subtype_by_map=TRUE)
celltype_map<-cell_activity_database$celltype_map

prefix<-outFile

if(!exists("obj")){
  obj<-read_object(parFile1)
}

resolutions=c(0.5, 1.0, 1.5)

curreduction=ifelse(by_harmony, "harmony", "pca")
obj <- FindNeighbors(object=obj, reduction=curreduction, dims=c(1:npcs), verbose=FALSE)
obj <- FindClusters(object=obj, verbose=FALSE, random.seed=random.seed, resolution=resolutions)

res_prefix<-paste0(assay, "_snn_res")
multi_res<-colnames(obj@meta.data)[grepl(paste0("^", res_prefix, ".+\\d$"), colnames(obj@meta.data))]

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
  
  g<-get_dim_plot_labelby(obj, label.by=cur_celltype, reduction="umap")
  g<-g+get_bubble_plot(obj, cur_res=NA, cur_celltype, bubblemap_file, assay="RNA", orderby_cluster = T)
  layout<-"
ABB
"
  g<-g+plot_layout(design = layout)
  png(paste0(prefix, ".", cur_celltype, ".png"), width=6600, height=2000, res=300)
  print(g)
  dev.off()
  
  g<-get_dim_plot(obj, group.by=cur_res, label.by=s_rawname, reduction="umap")
  g<-g+get_bubble_plot(obj, cur_res=cur_res, raw_celltype, bubblemap_file, assay="RNA", orderby_cluster = T)
  layout<-"
ABB
"
  g<-g+plot_layout(design = layout)
  png(paste0(prefix, ".", raw_celltype, ".seurat.png"), width=6600, height=2000, res=300)
  print(g)
  dev.off()
  
  g<-get_dim_plot(obj, group.by=cur_res, label.by=sname, reduction="umap")
  g<-g+get_bubble_plot(obj, cur_res=cur_res, cur_celltype, bubblemap_file, assay="RNA", orderby_cluster = T)
  layout<-"
ABB
"
  g<-g+plot_layout(design = layout)
  png(paste0(prefix, ".", cur_celltype, ".seurat.png"), width=6600, height=2000, res=300)
  print(g)
  dev.off()
  
  #draw_bubble_plot(obj, cur_res, cur_celltype, bubblemap_file, paste0(prefix, ".", cur_celltype, assay))
  
  save_highlight_cell_plot(paste0(prefix, ".", cur_celltype, ".cell.png"), obj, group.by = cur_celltype, reduction = "umap")

  if(file.exists(parSampleFile2)){
    new_layer = paste0(cur_celltype, "_summary")
    new_celltype<-as.character(obj@meta.data[,cur_celltype])
    in_new_celltype<-new_celltype %in% names(combined_ct)
    new_celltype[in_new_celltype]<-combined_ct[new_celltype[in_new_celltype]]
    obj[[new_layer]] = new_celltype

    g<-get_dim_plot_labelby(obj, label.by=new_layer, reduction="umap")
    g<-g+get_bubble_plot(obj, cur_res=NA, new_layer, bubblemap_file, assay="RNA", orderby_cluster = T)
    layout<-"
ABB
"
    g<-g+plot_layout(design = layout)
    png(paste0(prefix, ".", new_layer, ".png"), width=6600, height=2000, res=300)
    print(g)
    dev.off()
  }
}

saveRDS(obj@meta.data, paste0(prefix, ".meta.rds"))
write.csv(obj@meta.data, paste0(prefix, ".meta.csv"))

#saveRDS(obj, paste0(prefix, ".multires.rds"))
