rm(list=ls())
outFile='veteran'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/ravi_shah_projects/shengq2/20230327_rnaseq_veteran_hg38/md5_merge/result')

### Parameter setting end ###

md5_files = read.table("fileList1.txt", header=F, sep="\t")$V1
md5_values = lapply(md5_files, function(x){
  #if(file.exists(x)){
  return(readLines(x))
  #}
})
md5_values_all = unlist(md5_values)
writeLines(md5_values_all, paste0(outFile,".md5.txt"))
