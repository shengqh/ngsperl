source("split_samples_utils.r")
library(Seurat)
library(ggplot2)
library(scDemultiplex)
library(tictoc)

library(Seurat)
library(ggplot2)
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)
library(TailRank) # dbb # Beta-Binomial Distribution
library(edgeR)
library(dirmult)
library(MGLM) # ddirmn
library(parallel)
library("choisycutoff")

do_scDemultiplex<-function(fname, rdsfile, init_by_HTO_demux=0, cutoff_startval=0, valid_tags=NA) {
  print(fname, ":", rdsfile, " ...")

  output_prefix = paste0(fname, ".HTO")

  cur_tags = NULL
  if(!is.na(valid_tags)){
    cur_tags = strsplit(valid_tags, ",")
  }
  obj=read_hto(rdsfile, output_prefix, cur_tags)

  p.cut=0.001

  refine_rds<-paste0(output_prefix, ".scDemultiplex.refine.rds")
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
      print(paste0("starting ", fname, " by cutoff ..."))
      tic()
      obj<-demulti_cutoff(obj, output_prefix, cutoff_startval, mc.cores=nrow(obj))
      toc1=toc()
      init_column = "scDemultiplex_cutoff";
    }
    tic(paste0("refining ", fname, " ...\n"))
    obj<-demulti_refine(obj, p.cut, init_column=init_column, mc.cores=nrow(obj))
    toc2=toc()
    saveRDS(obj@meta.data, refine_rds)

    saveRDS(list("cutoff"=toc1, "refine"=toc2), paste0(output_prefix, ".scDemultiplex.tictoc.rds"))
  }else{
    obj@meta.data<-readRDS(refine_rds)
  }

  obj$HTO_classification<-obj$scDemultiplex
  obj$HTO_classification.global<-obj$scDemultiplex.global

  print(paste0("output result ..."))
  output_post_classification(obj, output_prefix)

  return(obj)
}

args = commandArgs(trailingOnly=TRUE)
print(paste0(args, collapse = " "))

if (length(args) == 0) {
  fname = "Adipose_9240"
  rdsfile = "/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples_preparation/result/Adipose_9240.hto.rds"
  init_by_HTO_demux = 0
  cutoff_startval = 0
  valid_tags = NA
}else{
  fname = args[1]
  rdsfile = args[2]
  init_by_HTO_demux = is_one(args[3])
  cutoff_startval = is_one(args[4])
  valid_tags = args[5]
}

obj = do_scDemultiplex( fname = fname, 
                        rdsfile = rdsfile, 
                        init_by_HTO_demux = init_by_HTO_demux, 
                        cutoff_startval = cutoff_startval, 
                        valid_tags = valid_tags)
