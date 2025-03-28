rm(list=ls()) 
sample_name='cvd_9a'
outFile='cvd_9a'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20241030_Kaushik_Amancherla_snRNAseq/20250211_T04_snRNA_hg38/raw_qc/result/cvd_9a')

### Parameter setting end ###

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RCurl)
library(knitr)
library(kableExtra)
library(data.table)

source("scRNA_func.r")
source("markerCode_filter.R")

options(future.globals.maxSize= 10779361280)

myoptions<-read_file_map("fileList2.txt", do_unlist=FALSE)
myoptions$mt_cutoff=as.numeric(myoptions$mt_cutoff)
myoptions$nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
myoptions$nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
myoptions$nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
myoptions$pca_dims=as.numeric(myoptions$pca_dims)
myoptions$Remove_MtRNA=is_one(myoptions$Remove_MtRNA)
myoptions$Remove_rRNA=is_one(myoptions$Remove_rRNA)

regress_by_percent_mt=is_one(myoptions$regress_by_percent_mt)
if(regress_by_percent_mt){
  vars.to.regress=c("percent.mt")
}else{
  vars.to.regress=NULL
}

Ensemblfile=myoptions$ensembl_gene_map_file
bubblemap_file = myoptions$bubblemap_file
bubblemap_width=to_numeric(myoptions$bubblemap_width, 4000)
bubblemap_height=to_numeric(myoptions$bubblemap_height, 2000)

Mtpattern=myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_Mt_rRNA=myoptions$Remove_MtRNA
resolution=as.numeric(myoptions$resolution)

by_sctransform=is_one(myoptions$by_sctransform)
use_sctransform_v2=is_one(myoptions$use_sctransform_v2)
output_object=is_one(myoptions$output_object)

species=myoptions$species # Hs or Mm

ignore_variable_genes=c()
if("ignore_variable_gene_file" %in% names(myoptions)){
  ignore_variable_gene_file = myoptions$ignore_variable_gene_file
  if(ignore_variable_gene_file != ""){
    ignore_variable_genes=readLines(ignore_variable_gene_file)
  }
}

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

Cutoffs=init_cutoffs(SampleInfos$SampleId, myoptions, parFile1)

miss_samples = SampleInfos$SampleId[!(SampleInfos$SampleId %in% Cutoffs$sample)]
if(length(miss_samples) > 0){
  stop(paste0("Samples ", paste0(miss_samples, collapse=","), " not in ", parFile1))
}

transpose=FALSE

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = species,
                             curated_markers_file = myoptions$curated_markers_file,
                             HLA_panglao5_file = myoptions$HLA_panglao5_file,
                             layer="Layer4",
                             remove_subtype_str = NULL,
                             combined_celltype_file = NULL)
tiers = ctdef$tiers
cell_activity_database<-ctdef$cell_activity_database
cellType = cell_activity_database$cellType

celltype_predictmethod="cta" # cta: cell activity; ora: over-represent

object.list<-list()

i=1
for (i in 1:nrow(SampleInfos)) {
  SampleInfo<-SampleInfos[i,]
  Cutoff<-Cutoffs[SampleInfo$SampleId,]
  info<-preprocess( SampleInfo = SampleInfo,
                    Cutoff = Cutoff,
                    cellType = cellType,
                    Mtpattern = Mtpattern,
                    rRNApattern = rRNApattern,
                    resolution = resolution, 
                    Remove_Mt_rRNA = Remove_Mt_rRNA,
                    celltype_predictmethod = celltype_predictmethod,
                    transpose = transpose,
                    hto_map = hto_map,
                    tag_tb = tag_tb,
                    Ensemblfile = Ensemblfile,
                    bubblemap_file = bubblemap_file,
                    bubblemap_width = bubblemap_width,
                    bubblemap_height = bubblemap_height,
                    species = species,
                    by_sctransform = by_sctransform,
                    use_sctransform_v2 = use_sctransform_v2,
                    output_object = output_object,
                    vars.to.regress = vars.to.regress,
                    ignore_variable_genes = ignore_variable_genes)
  
  object.list<-c(object.list, info)
}

saveRDS(object.list,file=paste0("objectlist.rds"))

stats<-lapply(object.list, function(x){unlist(x$preprocess)})
stats_df<-data.frame(do.call(rbind, stats))
colnames(stats_df)<-gsub("preprocess.","",colnames(stats_df))
write.table(stats_df, file="qc_filter_config.txt", sep="\t", row.names=F)

