rm(list=ls()) 
outFile='AG3669'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/cqs/softwares/ngsperl/lib/Pipeline/../scRNA/gliph2_config.txt'
parFile2='/data/h_gelbard_lab/projects/20220508_scRNA_3669/tcr_hla_data/result/AG3669.tcr.CD4.txt'
parFile3='/data/h_gelbard_lab/projects/20220508_scRNA_3669/tcr_hla_data/result/AG3669.tcr.CD8.txt'
parFile4='/data/h_gelbard_lab/projects/20220508_scRNA_3669/tcr_hla_data/result/AG3669.hla.txt'


setwd('/data/h_gelbard_lab/projects/20220508_scRNA_3669/tcr_hla_data_gliph2/result')

### Parameter setting end ###

library(data.table)
library(reshape2)

project_name <- outFile

config=read.table(parFile1, sep="=", stringsAsFactors = F, header=F)
rownames(config)<-config$V1

ref_files = read.table(parSampleFile1, sep="\t", stringsAsFactors = F, header=F)

tcr_files=list("CD4"=parFile2, "CD8"=parFile3)

tcr="CD4"
for(tcr in names(tcr_files)){
  tcr_file=tcr_files[tcr]
  tcr_ref = ref_files[ref_files$V3==tcr,]
  tcr_ref_map = unlist(split(tcr_ref$V1, tcr_ref$V2))
  config["out_prefix","V2"] = paste0(project_name, ".", tcr)
  config["cdr3_file","V2"] = tcr_file
  config["hla_file","V2"] = ifelse(exists("parFile4"), parFile4, "")
  config["refer_file","V2"] = tcr_ref_map['refer_file']
  config["v_usage_freq_file","V2"] = tcr_ref_map['v_usage_freq_file']
  config["cdr3_length_freq_file","V2"] = tcr_ref_map['cdr3_length_freq_file']
  config_file = paste0(project_name, ".", tcr, ".config.txt")
  write.table(config, config_file, sep="=", col.names = F, row.names=F, quote=F)

  cmd=paste0("irtools -c ", config_file)
  cat(cmd, "\n")
  system(cmd)
}
