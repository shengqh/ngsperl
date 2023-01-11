rm(list=ls()) 
outFile='P9112'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1='/data/h_gelbard_lab/projects/20221129_9112_scRNA_human/hto_samples.txt'
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20221129_9112_scRNA_human/hto_samples_HTODemux_all/result')

### Parameter setting end ###

source("split_samples_utils.r")
library(Seurat)
library(ggplot2)


files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

has_valid_tag = exists("parSampleFile4")
if(has_valid_tag){
  valid_tags = read.table(parSampleFile4, sep="\t", header=F)
  valid_tags_map = split(valid_tags$V1, valid_tags$V2)
}

idx=1
for(idx in c(1:length(files))){
  fname=names(files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  
  if(file.exists(output_file) & params$hto_ignore_exists){
    next
  }

  rdsfile=files[[idx]]
  cat(fname, ":", rdsfile, " ...\n")

  cur_tags = NULL
  if(has_valid_tag){
    if(fname %in% names(valid_tags_map)){
      cur_tags = unlist(valid_tags_map[fname])
    }else{
      warning(paste0(fname, " is not in valid tag map, check HTO_tags in configuration file, will continue to use all tags in this sample!"))
    }
  }
  obj=read_hto(rdsfile, output_prefix, cur_tags)

  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)

  obj$HTO_classification[obj$HTO_classification.global == "Doublet"] = "Doublet"

  output_post_classification(obj, output_prefix)
}
