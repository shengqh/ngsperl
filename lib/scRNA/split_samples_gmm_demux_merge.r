rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile5='fileList5.txt'
parFile1='/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples.txt'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/shengq2/20230115_combined_scRNA_hg38/hto_samples_gmm_demux/result')

### Parameter setting end ###

source("split_samples_utils.r")
library(Seurat)
library(ggplot2)

#devtools::install_github("shengqh/cutoff")
#install.packages("bbmle")
library(zoo)
library(reshape2)
library(gridExtra)
library(ggExtra)

full_files=read_map(parSampleFile1)
config_files=read_map(parSampleFile2)
hto_files=read_map(parSampleFile5)

params=read_map(parSampleFile3)
umap_min_dist=as.numeric(params$umap_min_dist)
umap_num_neighbors=as.numeric(params$umap_num_neighbors)

has_hto_samples = FALSE
if(exists('parSampleFile4')){
  has_hto_samples=TRUE
  hto_samples_tbl = read.table(parSampleFile4, sep="\t", header=F)
  hto_samples = split(hto_samples_tbl$V1, hto_samples_tbl$V2)
}
if(!has_hto_samples & (parFile1 != "")){
  has_hto_samples=TRUE
  hto_samples = read.table(parFile1, sep="\t", header=T)
}

idx=1
for(idx in c(1:length(full_files))){
  fname=names(full_files)[idx]
  output_prefix = paste0(fname, ".HTO")
  output_file=paste0(output_prefix, ".csv")
  
  #if(file.exists(output_file) & params$hto_ignore_exists){
  #  next
  #}

  full_file=full_files[[idx]]
  cat(fname, ":", full_file, " ...\n")

  config_file=config_files[[idx]]
  cat(fname, ":", config_file, " ...\n")

  hto_file=hto_files[[idx]]
  cat(fname, ":", hto_file, " ...\n")

  if (has_hto_samples){
    sample_tbl=subset(hto_samples, File==fname)
    cur_tags=sample_tbl$Tagname
    tag_sample_map=unlist(split(sample_tbl$Sample, sample_tbl$Tagname))
    tag_sample_map$Doublet<-"Doublet"
    tag_sample_map$Negative<-"Negative"
  }else{
    cur_tags = NULL
  }

  obj=read_hto(hto_file, output_prefix, cur_tags)
  
  tagnames=rownames(obj[["HTO"]])

  data <- FetchData(object=obj, vars=tagnames)
  write.csv(data, file=paste0(output_prefix, ".data.csv"))

  full <- read.csv(full_file, row.names=1, stringsAsFactor=F) 
  config <- read.csv(config_file, sep=',', header=F)
  config$V2 <- gsub('\\s+','', config$V2)

  negative_id<-config$V1[config$V2 == "negative"]
  singlet_clusters<-config[config$V2 %in% tagnames,]
  singlet_map=split(singlet_clusters$V2, singlet_clusters$V1)
  full$HTO_classification=unlist(lapply(full$Cluster_id, function(x){
    if (x == negative_id) {
      return("Negative")
    }
    if (as.character(x) %in% names(singlet_map)) {
      return(singlet_map[[as.character(x)]])
    }
    return("Doublet")
  }))
  
  full$HTO_classification.global=unlist(lapply(full$HTO_classification, function(x){
    if (!(x %in% c("Negative", "Doublet"))){
      return("Singlet")
    }else{
      return(x)
    }
  }))
  
  obj[["HTO_classification"]] = full$HTO_classification
  obj[["HTO_classification.global"]] = full$HTO_classification.global

  if(has_hto_samples){
    obj$orig.ident<-unlist(tag_sample_map[obj$HTO_classification])
  }

  saveRDS(obj, file=paste0(output_prefix, ".umap.rds"))
  
  output_post_classification(obj, output_prefix, umap_min_dist=umap_min_dist, umap_num_neighbors=umap_num_neighbors, tagnames=tagnames)
}
