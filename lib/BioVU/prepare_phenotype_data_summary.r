rm(list=ls()) 
outFile='phewas'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/celly_wanjalla_projects/20241209_phenotype_agd163k/phenotype_summary/result')

### Parameter setting end ###

library(data.table)
library(dplyr)

files=fread(parSampleFile1, header=FALSE, data.table=FALSE) |>
  dplyr::rename(file=1, phenotype=2)

res_df=data.frame()
idx=1
for(idx in c(1:nrow(files))){
  file<-files$file[idx]
  phenotype<-files$phenotype[idx]

  data <- fread(file, data.table=FALSE)
  data_res=data.frame("phenotype"=phenotype, "positive"=sum(data$Phenotype==1), "control"=sum(data$Phenotype==0))
  
  res_df=rbind(res_df, data_res)
}

write.csv(res_df, paste0(outFile, ".phenotype_summaey.csv"), row.names=FALSE)
