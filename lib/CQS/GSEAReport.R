rm(list=ls()) 
outFile='AG_integrated'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20220907_8566_project/seurat_sct_harmony_multires_03_choose_edgeR_inCluster_byCell_GSEA_Hs_report/result')

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

source('Pipeline.R')

files<-read.table(parSampleFile1, header=FALSE, as.is=TRUE)
if(nrow(files) == 1 & files$V2[1] == outFile){
  #single cell GSEA result
  files<-read.csv(files$V1[1])
  if(all(grepl("^_", files$compName))){
    files$compName<-gsub("_GSEA.rnk.*", "", basename(dirname(files$Folder)))
  }
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

