rm(list=ls()) 
outFile='GPA'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20230806_gpa_scRNA_hg38_merge_files/result')

### Parameter setting end ###

source("scRNA_func.r")

sample_map = read_file_map(parSampleFile1, do_unlist=FALSE)
h5_map = read_file_map(parSampleFile2)

merge_files<-function(sample_map, file_map, suffix){
  sample = names(sample_map)[1]
  for(sample in names(sample_map)){
    snames = sample_map[[sample]]
    cat(sample, ":", paste0(snames, collapse=", "), "\n")

    if(!all(snames %in% names(h5_map))){
      stop(paste0("Samples ", paste0(snames[!(snames %in% names(h5_map))], collapse=","), " not in ", parSampleFile2))
    }

    raw_objs = list()
    for(sname in snames){
      sfile = file_map[[sname]]
      cat("  read", sfile, "\n")
      counts = read_scrna_data(sfile)$counts
      colnames(counts)<-paste0(sname, "_", colnames(counts))
      sobj = CreateSeuratObject(counts = counts, project = sample)
      raw_objs[sname] = sobj
    }

    cat("merging all sub samples ... \n")
    rawobj <- merge(raw_objs[[1]], y = unlist(raw_objs[2:length(raw_objs)]), project = "sample")
    saveRDS(rawobj$RNA@counts, paste0(sample, suffix))
  }
}

merge_files(sample_map, h5_map, ".counts.filtered.rds")

if(parSampleFile3 != ''){
  raw_h5_map = read_file_map(parSampleFile3, do_unlist=FALSE)
  merge_files(sample_map, raw_h5_map, ".counts.raw.rds")
}
