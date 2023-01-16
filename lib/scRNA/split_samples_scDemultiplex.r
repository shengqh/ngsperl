rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1='/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples.txt'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples_scDemultiplex/result')

### Parameter setting end ###

source("split_samples_utils.r")
library(Seurat)
library(ggplot2)
library(scDemultiplex)
library(tictoc)

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

init_by_HTO_demux=is_one(params$init_by_HTO_demux)
cutoff_startval=ifelse(is.null(params$cutoff_startval), 0, ifelse(params$cutoff_startval == "", 0, as.numeric(params$cutoff_startval)))

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

  p.cut=0.001

  refine_rds<-paste0(fname, ".scDemultiplex.refine.rds")
  if(!file.exists(refine_rds)){

    if(init_by_HTO_demux){
      print(paste0("starting ", fname, " by HTODemux ..."))
      tic()
      obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
      toc1=toc()
      obj$HTO_classification[obj$HTO_classification.global == "Doublet"] = "Doublet"
      obj$HTODemux = obj$HTO_classification
      obj$HTODemux.global = obj$HTO_classification.global
      obj$HTO_classification = NULL
      obj$HTO_classification.global = NULL
      init_column = "HTODemux";
    }else{
      cutoff_rds<-paste0(fname, ".scDemultiplex.cutoff.startval_", cutoff_startval, ".rds")
      if(!file.exists(cutoff_rds)){
        print(paste0("starting ", fname, " by cutoff ..."))
        tic()
        obj<-demulti_cutoff(obj, output_prefix, cutoff_startval)
        toc1=toc()
        saveRDS(obj, cutoff_rds)
      }else{
        obj<-readRDS(cutoff_rds)
        toc1<-NA
      }
      init_column = "scDemultiplex_cutoff";
    }
    tic(paste0("refining ", fname, " ...\n"))
    obj<-demulti_refine(obj, p.cut, init_column=init_column)
    toc2=toc()
    saveRDS(obj@meta.data, refine_rds)

    saveRDS(list("cutoff"=toc1, "refine"=toc2), paste0(output_prefix, ".scDemultiplex.tictoc.rds"))
  }else{
    obj@meta.data<-readRDS(refine_rds)
  }

  obj$HTO_classification<-obj$scDemultiplex
  obj$HTO_classification.global<-obj$scDemultiplex.global

  cat(paste0("output result...\n"))
  output_post_classification(obj, output_prefix)
}
