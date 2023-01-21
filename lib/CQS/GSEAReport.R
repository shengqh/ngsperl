rm(list=ls()) 
outFile='crs'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_turner_lab/shengq2/20221206_7114_8822_scRNA_hg38/seurat_sct_merge_dr0.5_03_choose_edgeR_inCluster_bySample_GSEA_report/result')

### Parameter setting end ###

library(knitr)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(DT)
library(RCurl)
library(htmltools)
library(knitr)
library(kableExtra)

source('Functions.R')
source('Pipeline.R')

files<-read.table(parSampleFile1, header=FALSE, as.is=TRUE)
if(nrow(files) == 1 & files$V2[1] == outFile){
  #single cell GSEA result
  files<-read.csv(files$V1[1])
}else if ("task_name" %in% files$V2) {
  params_map<-split(files$V1, files$V2)
  task_name<-params_map$task_name
  gsea_file<-paste0(task_name, ".gsea.files.csv")
  files<-read.csv(gsea_file)
  files$Comparisons<-gsub(".gsea$", "", basename(dirname(files$Folder)))
}else{
  rownames(files)<-files$V2
}

vfiles = display_gsea(files, "", print_rmd = FALSE)
write.csv(vfiles, paste0("gsea_files.csv"))

