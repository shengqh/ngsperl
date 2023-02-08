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

library(Seurat)
library(ggplot2)
library(scDemultiplex)
library(tictoc)

source("split_samples_utils.r")

files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

init_by_HTODemux=is_one(params$init_by_HTODemux)
cutoff_startval=ifelse(is.null(params$cutoff_startval), 0.5, ifelse(params$cutoff_startval == "", 0.5, as.numeric(params$cutoff_startval)))

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

  cur_tags = NULL
  if(has_valid_tag){
    if(fname %in% names(valid_tags_map)){
      cur_tags = unlist(valid_tags_map[fname])
    }else{
      warning(paste0(fname, " is not in valid tag map, check HTO_tags in configuration file, will continue to use all tags in this sample!"))
    }
  }

  obj = do_scDemultiplex( fname=fname, 
                          rdsfile=rdsfile, 
                          output_prefix=output_prefix, 
                          init_by_HTODemux=init_by_HTODemux, 
                          cutoff_startval=cutoff_startval,
                          cur_tags=cur_tags)
                          
  saveRDS(obj, paste0(output_prefix, ".scDemultiplex.final.rds"))
}
