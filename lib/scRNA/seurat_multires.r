rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge/result/AK6383.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_multires/result')

### Parameter setting end ###

source("scRNA_func.r")
library(dplyr)
library(Seurat)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(rmdformats)
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
remove_subtype<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")
by_harmony<-reduction=="harmony"

bubblemap_file=myoptions$bubblemap_file

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

tiers<-read.table(myoptions$HLA_panglao5_file, sep="\t", header=T)

remove_subtype_of=remove_subtype
cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file, curated_markers_file=myoptions$curated_markers_file)

prefix<-outFile

if(!exists("obj")){
  obj<-read_object(parFile1)
}

resolutions=c(0.5, 1.0, 1.5)

curreduction=ifelse(by_harmony, "harmony", "pca")
obj <- FindNeighbors(object=obj, reduction=curreduction, dims=c(1:npcs), verbose=FALSE)
obj <- FindClusters(object=obj, verbose=FALSE, random.seed=random.seed, resolution=resolutions)

multi_res<-colnames(obj@meta.data)[grepl("_snn_res", colnames(obj@meta.data))]

cur_res = multi_res[1]
for(cur_res in multi_res){
  cat("  ", cur_res, "\n")
  data.norm=get_seurat_average_expression(obj, cur_res)
  
  predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
  
  new.cluster.ids<-names(predict_celltype$max_cta)
  names(new.cluster.ids) <- colnames(data.norm)
  
  cur_celltype<-paste0(cur_res, "_celltype")
  res_values=unlist(FetchData(obj, cur_res))
  obj[[cur_celltype]] = unlist(new.cluster.ids[res_values])
}

multi_cts<-colnames(obj@meta.data)[grepl("_celltype", colnames(obj@meta.data))]

res_df<-data.frame("resolution"=resolutions, "cluster"=multi_res, "celltype"=multi_cts)
write.csv(res_df, file=paste0(outFile, ".resolutions.csv"), row.names=F)

#umaplist<-RunMultipleUMAP(obj, curreduction=curreduction, cur_pca_dims=c(1:npcs))

umaplist<-RunMultipleUMAP(obj, nn=c(30,20,10), min.dist=c(0.3,0.3,0.3), curreduction=curreduction, cur_pca_dims=c(1:npcs))
obj<-umaplist$obj
umap_names<-umaplist$umap_names
rm(umaplist)

meta<-obj@meta.data

cur_celltype<-multi_cts[1]
for(cur_celltype in multi_cts){
  cat(cur_celltype, "\n")
  
  cur_res = gsub("_celltype", "", cur_celltype)
  meta<-meta[order(meta[,cur_res]),]  
  seurat_celltype<-paste0(meta[,cur_res], ":", meta[,cur_celltype])
  seurat_celltype<-factor(seurat_celltype, levels=unique(seurat_celltype))
  sname = paste0("seurat_", cur_celltype)
  meta[,sname] = seurat_celltype
  meta<-meta[colnames(obj),]
  obj@meta.data = meta
  
  g<-NULL
  for(umap_name in umap_names){
    g1<-get_dim_plot_labelby(obj, label.by=cur_celltype, reduction=umap_name) + ggtitle(umap_name)
    #g1<-DimPlot(obj, group.by = cur_celltype, reduction=umap_name, label=T, repel=T) + ggtitle(umap_name) + guides(fill=guide_legend(ncol=1))
    if(is.null(g)){
      g<-g1
    }else{
      g<-g+g1
    }
  }
  g<-g+get_bubble_plot(obj, cur_res=NA, cur_celltype, bubblemap_file, assay="RNA")
  layout<-"
ABC
DDD
"
  g<-g+plot_layout(design = layout)
  png(paste0(prefix, ".", cur_celltype, ".png"), width=as.numeric(myoptions$plot_width), height=as.numeric(myoptions$plot_height), res=300)
  print(g)
  dev.off()

  g<-NULL
  for(umap_name in umap_names){
    g2<-get_dim_plot(obj, group.by=cur_res, label.by=sname, reduction=umap_name) + ggtitle(umap_name)
    #g2<-DimPlot(obj, group.by = sname, reduction=umap_name, label=T, repel=T) + guides(color=guide_legend(ncol=1))
    if(is.null(g)){
      g<-g2
    }else{
      g<-g+g2
    }
  }

  g<-g+get_bubble_plot(obj, cur_res=cur_res, cur_celltype, bubblemap_file, assay="RNA")
  layout<-"
ABC
DDD
"
  g<-g+plot_layout(design = layout)
  png(paste0(prefix, ".", cur_celltype, ".seurat.png"), width=as.numeric(myoptions$plot_width), height=as.numeric(myoptions$plot_height), res=300)
  print(g)
  dev.off()
  
  #draw_bubble_plot(obj, cur_res, cur_celltype, bubblemap_file, paste0(prefix, ".", cur_celltype, assay))
  
  umap_name<-umap_names[1]
  
  save_highlight_cell_plot(paste0(prefix, ".", cur_celltype, ".cell.png"), obj, group.by = cur_celltype, reduction = umap_name)
}

saveRDS(obj@meta.data, paste0(prefix, ".meta.rds"))
write.csv(obj@meta.data, paste0(prefix, ".meta.csv"))

#saveRDS(obj, paste0(prefix, ".multires.rds"))
