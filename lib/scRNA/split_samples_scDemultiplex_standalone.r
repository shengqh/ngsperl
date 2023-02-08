library(Seurat)
library(ggplot2)
library(scDemultiplex)
library(tictoc)
library(stringr)

source("split_samples_utils.r")

args = commandArgs(trailingOnly=TRUE)
print(paste0(args, collapse = " "))

if (length(args) == 0) {
  fname = "Adipose_9240"
  rdsfile = "/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples_preparation/result/Adipose_9240.hto.rds"
  output_prefix = "Adipose_9240.HTO"
  init_by_HTO_demux = 0
  cutoff_startval = 0.5
  cur_tags = NULL
}else{
  fname = args[1]
  rdsfile = args[2]
  output_prefix = args[3]
  init_by_HTO_demux = is_one(args[4])
  cutoff_startval = ifelse(length(args) >= 5, as.numeric(args[5]), 0.5)
  cur_tags = ifelse(length(args) >= 6, str_split(args[6], ','), NULL)
}

obj = do_scDemultiplex( fname=fname, 
                        rdsfile=rdsfile, 
                        output_prefix=output_prefix, 
                        init_by_HTODemux=init_by_HTODemux, 
                        cutoff_startval=cutoff_startval,
                        cur_tags=cur_tags)

saveRDS(obj, paste0(output_prefix, ".scDemultiplex.final.rds"))

