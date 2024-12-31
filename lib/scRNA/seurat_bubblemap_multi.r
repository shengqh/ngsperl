rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.final.rds'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose_bubblemap/result')

### Parameter setting end ###

source("scRNA_func.r")
library("Seurat")
library("readxl")
library("ggplot2")

bubble_files=read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
cluster_name = ifelse(is.null(myoptions$cluster_name), "seurat_clusters", ifelse(myoptions$cluster_name=="", "seurat_clusters", myoptions$cluster_name))
celltype_name = ifelse(is.null(myoptions$celltype_name), "cell_type", ifelse(myoptions$celltype_name=="", "cell_type", myoptions$celltype_name))
summary_layer = ifelse(is.null(myoptions$summary_layer), "layer4", ifelse(myoptions$summary_layer=="", "layer4", myoptions$summary_layer))

has_summary=summary_layer != celltype_name

if(!exists('obj')){
  obj<-read_object(parFile1, parFile2)
}

cur_folder = getwd()
tmp_folder = paste0(cur_folder, "/details")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}

draw_figure<-function(subobj, group1, group2, bubblemap_file, png_prefix, rotate.title, by_group1, width, min_height, height_per_entry, height_additional_space){
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

  if(group2 == ""){
    height=get_dot_height_num(length(unique(subobj$dump_cluster_c)), min_height, height_per_entry, height_additional_space)
    ggsave( paste0(png_prefix, ".png"), 
            g,
            width=width, 
            height=min(50000, get_dot_height_num(length(unique(subobj$dump_cluster_c)), min_height, height_per_entry, height_additional_space)),
            units="px",
            dpi=300,
            bg="white",
            limitsize=FALSE)
  }

  if(by_group1){
    galldata = g$data
    cts=unique(unlist(subobj[[group1]]))
    ct<-cts[1]
    for(ct in cts){
      c_meta<-meta[meta[,group1]==ct,]
      g_data<-galldata[galldata$id %in% c_meta$dump_cluster_c,]
      g$data<-g_data
      ggsave( paste0(png_prefix, ".", celltype_to_filename(ct), ".png"), 
              g,
              width=width, 
              height=min(50000, get_dot_height_num(length(unique(g_data$id)), min_height, height_per_entry, height_additional_space)),
              units="px",
              dpi=300,
              bg="white",
              limitsize=FALSE)
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
    cat("cell_type_pattern: ", cell_type_pattern, "\n")
    parts=unlist(strsplit(cell_type_pattern, ":"))
    if(length(parts) == 1){
      pattern_col=celltype_name
      pattern_val=parts[1]
    }else{
      pattern_col=parts[1]
      pattern_val=parts[2]
    }
    if(!pattern_col %in% colnames(obj@meta.data)){
      stop(paste0("Column ", pattern_col, " not found in meta data"))
    }
    cells = colnames(obj)[grepl(pattern_val, unlist(obj[[pattern_col]]), ignore.case = ignore.case)]
    subobj=subset(obj, cells=cells)
  }else{
    subobj=obj
  }

  draw_figure(subobj = subobj, 
              group1 = celltype_name, 
              group2 = "", 
              bubblemap_file = bubblemap_file, 
              png_prefix = paste0(tmp_folder, "/", outFile, ".", bn, ".bubblemap.ct"), 
              rotate.title = rotate.title, 
              by_group1 = FALSE, 
              width = width, 
              min_height = min_height, 
              height_per_entry = height_per_entry, 
              height_additional_space = height_additional_space)

  if(length(unique(subobj@meta.data[,celltype_name])) != length(unique(subobj@meta.data[,cluster_name]))){
    draw_figure(subobj = subobj, 
                group1 = celltype_name, 
                group2 = cluster_name, 
                bubblemap_file = bubblemap_file, 
                png_prefix = paste0(tmp_folder, "/", outFile, ".", bn, ".bubblemap.ct_cluster"), 
                rotate.title = rotate.title, 
                by_group1 = TRUE, 
                width = width, 
                min_height = min_height, 
                height_per_entry = height_per_entry, 
                height_additional_space = height_additional_space)
  }

  draw_figure(subobj = subobj, 
              group1 = celltype_name, 
              group2 = "orig.ident", 
              bubblemap_file = bubblemap_file, 
              png_prefix = paste0(tmp_folder, "/", outFile, ".", bn, ".bubblemap.ct_ident"), 
              rotate.title = rotate.title, 
              by_group1 = TRUE, 
              width = width, 
              min_height = min_height, 
              height_per_entry = height_per_entry, 
              height_additional_space = height_additional_space)

  draw_figure(subobj = subobj, 
              group1 = "orig.ident", 
              group2 = celltype_name, 
              bubblemap_file = bubblemap_file, 
              png_prefix = paste0(tmp_folder, "/", outFile, ".", bn, ".bubblemap.ident_ct"), 
              rotate.title = rotate.title, 
              by_group1 = TRUE, 
              width = width, 
              min_height = min_height, 
              height_per_entry = height_per_entry, 
              height_additional_space = height_additional_space)
}
