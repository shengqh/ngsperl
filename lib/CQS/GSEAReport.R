rm(list=ls()) 
outFile='SADIE_adipose'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240304_mona_scRNA_SADIE/20240305_DE_fold1.2_adipose/files_edgeR_inCluster_bySample_GSEA_Hs_test_report/result')

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
  files$GseaCategory = basename(files$GseaCategory)
}else if ("task_name" %in% files$V2) {
  params_map<-split(files$V1, files$V2)
  task_name<-params_map$task_name
  gsea_file<-paste0(task_name, ".gsea.files.csv")
  files<-read.csv(gsea_file)
  files$Comparisons<-gsub(".gsea$", "", basename(dirname(files$Folder)))
}else if(all(grepl("gsea.files.csv$", files$V1))){
  flist=files
  files<-NULL
  i=1
  for(i in 1:nrow(flist)){
    f<-read.csv(flist$V1[i])
    f$compName=basename(f$compName)
    f$Comparisons<-gsub("_GSEA.rnk.*", "", basename(dirname(dirname(f$Folder))))
    files<-rbind(files, f)
  }
}else{
  rownames(files)<-files$V2
}

gsea_folder = paste0(outFile, ".gsea/")
dir.create(gsea_folder, showWarnings = FALSE)

vfiles = display_gsea(files, gsea_folder, print_rmd = FALSE)
write.csv(vfiles, paste0("gsea_files.csv"))
