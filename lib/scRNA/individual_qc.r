rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230115_combined_scRNA_hg38/individual_qc/result')

### Parameter setting end ###

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RCurl)
library(knitr)
library(kableExtra)

source("scRNA_func.r")
source("markerCode_filter.R")

option_table<-read.table("fileList2.txt", sep="\t", stringsAsFactors = F)
myoptions<-split(option_table$V1, option_table$V2)
myoptions$mt_cutoff=as.numeric(myoptions$mt_cutoff)
myoptions$nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
myoptions$nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
myoptions$nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
myoptions$pca_dims=as.numeric(myoptions$pca_dims)
myoptions$resolution=as.numeric(myoptions$resolution)
myoptions$Remove_MtRNA=is_one(myoptions$Remove_MtRNA)
myoptions$Remove_rRNA=is_one(myoptions$Remove_rRNA)
Mtpattern=myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_Mt_rRNA=myoptions$Remove_MtRNA
resolution=as.numeric(myoptions$resolution)

species=myoptions$species # Hs or Mm

#raw data file locations
SampleInfos<-read.table("fileList1.txt", stringsAsFactors = F)
colnames(SampleInfos)<-c("countfile", "SampleId")

if(file.exists('fileList3.txt')){
  hto_samples<-read.table('fileList3.txt', stringsAsFactors = F)
  hto_map=split(hto_samples$V1, hto_samples$V2)
  
  tag_tb<-read.table(myoptions$hto_sample_file, sep="\t", stringsAsFactors = F, header=T)
}else{
  hto_map=list()
  tag_tb=NULL
}

##cutoff dataframe for each sample to filter empty droplets
Cutoffs<-data.frame(nFeature_cutoff_min=myoptions$nFeature_cutoff_min ,
                    nFeature_cutoff_max=myoptions$nFeature_cutoff_max,
                    nCount_cutoff=myoptions$nCount_cutoff, 
                    mt_cutoff=myoptions$mt_cutoff, 
                    cluster_remove=c(""),stringsAsFactors = F)

#every sample share the same cutoff, otherwise put every parameter in the Cutoff dataframe
if (dim(Cutoffs)[1]==1){
  Cutoffs<-data.frame(lapply(Cutoffs, rep, nrow(SampleInfos)),stringsAsFactors = F)
}

transpose=FALSE

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = species,
                             curated_markers_file = myoptions$curated_markers_file,
                             HLA_panglao5_file = myoptions$HLA_panglao5_file,
                             layer="Layer4",
                             remove_subtype_str = "",
                             combined_celltype_file = NULL)
tiers = ctdef$tiers
cell_activity_database<-ctdef$cell_activity_database

celltype_predictmethod="cta" # cta: cell activity; ora: over-represent

object.list<-list()

i=1
for (i in 1:nrow(SampleInfos)) {
  SampleInfo<-SampleInfos[i,]
  Cutoff<-Cutoffs[i,]
  bubble_file=myoptions$bubblemap_file
  Ensemblfile=NULL
  info<-preprocess( SampleInfo = SampleInfo,
                    Cutoff = Cutoff,
                    cellType = cell_activity_database$cellType,
                    Mtpattern = Mtpattern,
                    rRNApattern = rRNApattern,
                    resolution = resolution, 
                    Remove_Mt_rRNA = Remove_Mt_rRNA,
                    celltype_predictmethod = celltype_predictmethod,
                    transpose = transpose,
                    hto_map = hto_map,
                    tag_tb = tag_tb,
                    Ensemblfile = Ensemblfile,
                    bubble_file = bubble_file)
  
  object.list<-c(object.list, info)
}

saveRDS(object.list,file=paste0("objectlist.rds"))

stats<-lapply(object.list, function(x){unlist(x$preprocess)})
stats_df<-data.frame(do.call(rbind, stats))
colnames(stats_df)<-gsub("preprocess.","",colnames(stats_df))
write.table(stats_df, file="qc_filter_config.txt", sep="\t", row.names=F)
