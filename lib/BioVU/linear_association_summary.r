rm(list=ls()) 
outFile='CFTR_rs113993960'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/h_cqs/shengq2/biovu/phecode_data/phecode_definitions1.2.csv'
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/jennifer_pietenpol_projects/20241011_PheWAS_of_immune_checkpoint_genes/20241107_CFTR_rs113993960/linear_association_summary/result')

### Parameter setting end ###

library(data.table)
library(dplyr)

files=fread(parSampleFile1, data.table=FALSE) |>
  dplyr::rename(file=1, phecode=2)

res_df=data.frame()
idx=1
for(idx in c(1:nrow(files))){
  file<-files$file[idx]
  phecode<-files$phecode[idx]

  #Load
  data <- fread(file, data.table=FALSE)
  data_res=data[2,,drop=FALSE]    
  data_res$V1[1]=phecode

  res_df=rbind(res_df, data_res)
}

res_df=res_df |>
  dplyr::rename(phecode=V1)

phecodes_def=fread(parFile1, data.table=FALSE) |>
  dplyr::select(phecode, phenotype) |>
  dplyr::distinct()

final_df=merge(res_df, phecodes_def, by="phecode", all.x=TRUE) |>
  dplyr::select(phecode, phenotype, everything()) |>
  dplyr::rename(pvalue=6)
final_df$padj=p.adjust(final_df$pvalue, method="fdr")
final_df = final_df[order(final_df$pvalue),]

write.csv(final_df, paste0(outFile, ".linear_association.csv"), row.names=FALSE)

