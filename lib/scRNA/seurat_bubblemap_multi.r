rm(list=ls()) 
outFile='P9551'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/jbrown_lab/projects/20230329_9551_flynn_hg38_liver_seurat/seurat_merge/result/P9551.final.rds'
parFile2='/nobackup/jbrown_lab/projects/20230329_9551_flynn_hg38_liver_seurat/seurat_merge_dr0.5_01_call/result/P9551.scDynamic.meta.rds'
parFile3=''


setwd('/nobackup/jbrown_lab/projects/20230329_9551_flynn_hg38_liver_seurat/seurat_merge_dr0.5_01_call_bubblemap_final/result')

### Parameter setting end ###

source("scRNA_func.r")
library("Seurat")
library("readxl")
library("ggplot2")

bubble_files=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
cluster_name = ifelse(is.null(myoptions$cluster_name) | myoptions$cluster_name=="", "seurat_clusters", myoptions$cluster_name)
celltype_name = ifelse(is.null(myoptions$celltype_name) | myoptions$celltype_name=="", "cell_type", myoptions$celltype_name)

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2)
}

cur_folder = getwd()
tmp_folder = paste0(cur_folder, "/details")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}

draw_figure<-function(subobj, group1, group2, bubblemap_file, png_prefix, rotate.title, tmp_folder, by_group1, width, min_height, height_per_entry, height_additional_space){
  if(group2 == ""){
    subobj$dump_cluster = unlist(subobj[[group1]])
  }else{
    subobj$dump_cluster = paste0(unlist(subobj[[group1]]), ":", unlist(subobj[[group2]]))
  }

  dc_tbl<-table(subobj$dump_cluster)
  dc_tbl_c<-paste0(names(dc_tbl), "(", dc_tbl, ")")
  names(dc_tbl_c)<-names(dc_tbl)
  subobj$dump_cluster_c<-dc_tbl_c[subobj$dump_cluster]

  meta<-subobj@meta.data
  if(group2 == ""){
    meta<-meta[order(meta[[group1]]), ]
  }else{
    meta<-meta[order(meta[[group1]], meta[[group2]]), ]
  }

  subobj$dump_cluster_c<-factor(subobj$dump_cluster_c, unique(meta$dump_cluster_c))

  g=get_bubble_plot(subobj, 
    NULL, 
    NULL, 
    bubblemap_file, 
    assay="RNA", 
    orderby_cluster=TRUE, 
    rotate.title=rotate.title, 
    group.by="dump_cluster_c",
    species=myoptions$species) + scale_color_gradient2(low="blue", mid="yellow", high="red")


  png(paste0(png_prefix, ".png"), width=width, height=get_dot_height_num(length(unique(subobj$dump_cluster_c)), min_height, height_per_entry, height_additional_space),res=300)
  print(g)
  dev.off()

  if(by_group1){
    galldata = g$data
    cts=unique(unlist(subobj[[group1]]))
    ct<-cts[1]
    for(ct in cts){
      c_meta<-meta[meta[,group1]==ct,]
      g_data<-galldata[galldata$id %in% c_meta$dump_cluster_c,]
      g$data<-g_data
      png(paste0(tmp_folder, "/", png_prefix, ".", celltype_to_filename(ct), ".png"), width=width, height=get_dot_height_num(length(unique(g_data$id)), min_height, height_per_entry, height_additional_space),res=300)
      print(g)
      dev.off()
    }
  }
}

samples=unique(as.character(obj$orig.ident))
samples=samples[order(samples)]
writeLines(samples, paste0(outFile, ".sample.txt"))

celltypes=names(table(obj[[celltype_name]]))
writeLines(celltypes, paste0(outFile, ".cell_type.txt"))

bnames=unique(bubble_files$V3)
bn=bnames[1]
for(bn in bnames){
  cat(bn, "\n")

  bn_params<-bubble_files[bubble_files$V3==bn,]
  bn_map<-split(bn_params$V1, bn_params$V2)

  cell_type_pattern=bn_map$cell_type_pattern
  bubblemap_file=bn_map$file
  ignore.case=is_one(bn_map$ignore.case)
  rotate.title= is_one(bn_map$rotate.title)
  width=as.numeric(bn_map$width)
  min_height=to_numeric(bn_map$min_height, 1500)
  height_per_entry=to_numeric(bn_map$height_per_entry, 80)
  height_additional_space=to_numeric(bn_map$height_additional_space, 800)

  #height=as.numeric(bn_map$height)
  
  if(is.null(cell_type_pattern)){
    subobj=obj
  }else if (cell_type_pattern != "" & cell_type_pattern != "*"){
    cells = colnames(obj)[grepl(cell_type_pattern, unlist(obj[[celltype_name]]), ignore.case = ignore.case)]
    subobj=subset(obj, cells=cells)
  }else{
    subobj=obj
  }

  draw_figure(subobj, celltype_name, "", bubblemap_file, paste0(outFile, ".", bn, ".bubblemap.ct"), rotate.title, tmp_folder, by_group1 = FALSE, width, min_height, height_per_entry, height_additional_space)

  draw_figure(subobj, celltype_name, cluster_name, bubblemap_file, paste0(outFile, ".", bn, ".bubblemap.ct_cluster"), rotate.title, tmp_folder, TRUE, width, min_height, height_per_entry, height_additional_space)

  draw_figure(subobj, celltype_name, "orig.ident", bubblemap_file, paste0(outFile, ".", bn, ".bubblemap.ct_ident"), rotate.title, tmp_folder, TRUE, width, min_height, height_per_entry, height_additional_space)

  draw_figure(subobj, "orig.ident", celltype_name, bubblemap_file, paste0(outFile, ".", bn, ".bubblemap.ident_ct"), rotate.title, tmp_folder, TRUE, width, min_height, height_per_entry, height_additional_space)
}
