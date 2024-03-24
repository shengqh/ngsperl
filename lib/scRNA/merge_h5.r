rm(list=ls()) 
outFile='iSGS_cell_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240320_scRNA_iSGS_cell_atlas/cellbender_04_pool/result')

### Parameter setting end ###

source("scRNA_func.r")
library(DropletUtils)

sample_map = read_file_map(parSampleFile1, do_unlist=FALSE)
h5_map = read_file_map(parSampleFile2)
myoptions=read_file_map(parSampleFile4, do_unlist=FALSE)

merge_files<-function(sample_map, file_map, outFile, suffix){
  df = NULL
  sample = names(sample_map)[1]
  for(sample in names(sample_map)){
    snames = sample_map[[sample]]
    cat(sample, ":", paste0(snames, collapse=", "), "\n")

    if(!all(snames %in% names(file_map))){
      stop(paste0("Samples ", paste0(snames[!(snames %in% names(file_map))], collapse=","), " not in ", parSampleFile2))
    }

    raw_objs = list()
    sname=snames[1]
    for(sname in snames){
      sfile = file_map[[sname]]
      cat("  read", sfile, "\n")
      counts = read_scrna_data(sfile)$counts
      colnames(counts)<-paste0(sname, "_", colnames(counts))
      sobj = CreateSeuratObject(counts = counts, project = sample)
      raw_objs[sname] = sobj
    }

    cat("  merging all sub samples ... \n")
    rawobj <- merge(raw_objs[[1]], y = unlist(raw_objs[2:length(raw_objs)]), project = "sample")

    final_file=paste0(sample, suffix)
    cat("  saving to", final_file, "...\n")
    write10xCounts( final_file, rawobj$RNA@counts)   

    df=rbind(df, data.frame(sample=sample, ncells=ncol(rawobj$RNA@counts), ngenes=nrow(rawobj$RNA@counts))) 
  }
  write.table(df, paste0(outFile, suffix, ".txt"), row.names=FALSE, sep="\t", quote=FALSE)
}

merge_files(sample_map, h5_map, outFile, myoptions$suffix)

if(parSampleFile3 != ''){
  raw_h5_map = read_file_map(parSampleFile3, do_unlist=FALSE)
  merge_files(sample_map, raw_h5_map, outFile, gsub(".h5", ".raw.h5", myoptions$suffix))
}
