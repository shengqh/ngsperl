source("split_samples_utils.r")
library(Seurat)
library(ggplot2)


files_lines=read.table(parSampleFile1, sep="\t")
files=split(files_lines$V1, files_lines$V2)

cutoffs_lines=read.table(parSampleFile2, sep="\t")
cutoffs=split(cutoffs_lines$V1, cutoffs_lines$V2)

params_lines=read.table(parSampleFile3, sep="\t")
params=split(params_lines$V1, params_lines$V2)
params$hto_ignore_exists=ifelse(params$hto_ignore_exists=="0", FALSE, TRUE)

idx=3
for(idx in c(1:length(files))){
  fname=names(files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  
  if(file.exists(output_file) & params$hto_ignore_exists){
    next
  }

  h5file=files[[idx]]
  cat(fname, ":", h5file, " ...\n")

  obj=read_hto(h5file, output_prefix)

  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)

  obj$HTO_classification[obj$HTO_classification.global == "Doublet"] = "Doublet"

  output_post_classification(obj, output_prefix)
}
