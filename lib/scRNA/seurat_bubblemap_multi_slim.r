rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/seurat_sct2_merge/result/iSGS_cell_atlas.final.rds'
parFile2='/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.2_1_call_harmony/result/iSGS_cell_atlas.scDynamic.meta.rds'
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/seurat_sct2_merge_dr0.2_1_call_harmony_bubble_files/result')

### Parameter setting end ###

source("scRNA_func.r")
library("Seurat")
library("readxl")
library("ggplot2")

bubble_files=read_file_map(parSampleFile1)

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
cluster_name = ifelse(is.null(myoptions$cluster_name) | myoptions$cluster_name=="", "seurat_clusters", myoptions$cluster_name)
celltype_name = ifelse(is.null(myoptions$celltype_name) | myoptions$celltype_name=="", "cell_type", myoptions$celltype_name)
suffix = myoptions$suffix

file_dir=paste0(outFile, suffix)
dir.create(file_dir, showWarning=FALSE)
file_prefix=file.path(file_dir, outFile)

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2)
}

draw_figure<-function(subobj, group1, group2, bubblemap_file, png_file, rotate.title, tmp_folder, by_group1, width, min_height, height_per_entry, height_additional_space){
  if(group2 == ""){
    subobj@meta.data$dump_cluster = unlist(subobj[[group1]])
  }else{
    subobj@meta.data$dump_cluster = paste0(unlist(subobj[[group1]]), ":", unlist(subobj[[group2]]))
  }

  dc_tbl<-table(subobj$dump_cluster)
  dc_tbl_c<-paste0(names(dc_tbl), "(", dc_tbl, ")")
  names(dc_tbl_c)<-names(dc_tbl)
  subobj@meta.data$dump_cluster_c<-dc_tbl_c[subobj$dump_cluster]

  meta<-subobj@meta.data
  if(group2 == ""){
    meta<-meta[order(meta[[group1]]), ]
  }else{
    meta<-meta[order(meta[[group1]], meta[[group2]]), ]
  }

  subobj@meta.data$dump_cluster_c<-factor(subobj$dump_cluster_c, unique(meta$dump_cluster_c))

  g=get_bubble_plot(subobj, 
    NULL, 
    NULL, 
    bubblemap_file, 
    assay="RNA", 
    orderby_cluster=TRUE, 
    rotate.title=rotate.title, 
    group.by="dump_cluster_c",
    species=myoptions$species) + scale_color_gradient2(low="blue", mid="yellow", high="red")

  height=get_dot_height_num(length(unique(subobj$dump_cluster_c)), min_height, height_per_entry, height_additional_space)
  ggsave(png_file, g, width=width, height=height, unit="px", dpi=300, bg="white")
}

celltypes=names(table(obj[[celltype_name]]))
writeLines(celltypes, paste0(file_prefix, ".cell_type.txt"))

file_tbl=data.frame(data.frame(Name=names(bubble_files), ct_png="", ct_cluster_png=""), row.names=1)
bnames=names(bubble_files)
bn=bnames[1]
for(bn in bnames){
  cat(bn, "\n")

  bubblemap_file=bubble_files[[bn]]

  rotate.title= TRUE
  width=5000
  min_height=1500
  height_per_entry=80
  height_additional_space=800

  ct_png_file=paste0(file_prefix, ".bubblemap.", bn, ".ct.png")
  draw_figure(obj, celltype_name, "", bubblemap_file, ct_png_file, rotate.title, tmp_folder, by_group1 = FALSE, width, min_height, height_per_entry, height_additional_space)
  file_tbl[bn, "ct_png"]=ct_png_file

  ct_cluster_png_file=paste0(file_prefix, ".bubblemap.", bn, ".ct_cluster.png")
  draw_figure(obj, celltype_name, cluster_name, bubblemap_file, ct_cluster_png_file, rotate.title, tmp_folder, TRUE, width, min_height, height_per_entry, height_additional_space)
  file_tbl[bn, "ct_cluster_png"]=ct_cluster_png_file
}

write.csv(file_tbl, paste0(file_prefix, ".bubblemap.csv"), row.names=T, quote=F)
