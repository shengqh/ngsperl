rm(list=ls())
outFile='8857_human'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''

setwd('/nobackup/jbrown_lab/projects/20230404_atac_8857/encode_atacseq_local_croo_qc/result')

### Parameter setting end ###

library(data.table)

options(future.globals.maxSize= 10779361280)

option_tbl=fread(parSampleFile1, data.table=F, header=F)
myoptions = split(option_tbl$V1, option_tbl$V2)

outFile = myoptions$task_name

json_files = fread(parSampleFile2, data.table=F, header=F)
json_map = unlist(split(json_files$V1, json_files$V2))

sample_names = json_files$V2
sample_name = sample_names[1]
for(sample_name in sample_names){
  json_file=json_map[sample_name]
  #handle each json file
}

#generate table and figures for report, use outFile as prefix
