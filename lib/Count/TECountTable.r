rm(list=ls()) 
outFile='PP_6999_mouse_aorta_lmna'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/brown_lab/references/mm10_LMNA_Cas9/genes/genes.gtf.map'
parFile2='/nobackup/brown_lab/projects/2021/20211008_rnaseq_6999_PP_mouse_aorta_lmna/genetable/result/PP_6999_mouse_aorta_lmna.count'
parFile3=''


setwd('/nobackup/brown_lab/projects/2021/20211008_rnaseq_6999_PP_mouse_aorta_lmna/te_genetable/result')

### Parameter setting end ###

library(data.table)
library(dplyr)

fileList1 = fread(parSampleFile1, header=FALSE, data.table=FALSE) |>
  dplyr::rename("File"="V1", "Sample"="V2")

namemap = fread(parFile1, header=TRUE, data.table=FALSE) |>
  dplyr::select(gene_id, gene_name, gene_biotype)

final_df = NULL

i=1
for (i in 1:nrow(fileList1)){
  file = fileList1$File[i]
  sample = fileList1$Sample[i]

  cat("Processing ", file, " ... \n")
  df = fread(file, header=TRUE, data.table=FALSE)
  colnames(df)=c("gene_id", sample)

  if(is.null(final_df)){
    final_df = merge(df, namemap, by="gene_id", all.x=TRUE) |>
      dplyr::mutate(gene_biotype=ifelse(is.na(gene_name), "transposable elements", gene_biotype)) |>
      dplyr::select(gene_id, gene_name, gene_biotype, everything())

    final_df$gene_name[is.na(final_df$gene_name)] = final_df$gene_id[is.na(final_df$gene_name)]
  }else{
    final_df = merge(final_df, df, by="gene_id", all.x=TRUE)
  }
}

data_df=final_df |>
  dplyr::select(-gene_id, -gene_name, -gene_biotype)

has_reads=apply(data_df, 1, sum) > 0

valid_df=final_df[has_reads,] |>
  dplyr::rename("Feature"="gene_id", 
                "Feature_gene_name"="gene_name", 
                "Feature_gene_biotype"="gene_biotype") |>
  dplyr::select(Feature, Feature_gene_biotype, Feature_gene_name, everything())

if(parFile2 != ""){
  # For transposable elements, we will use the gene count from TEcount output which allows multimapping reads in STAR
  # But for gene, we will use the gene count from normal STAR/featureCount output which allows unique reads only
  cat("Extracting gene count from ", parFile2, " ... \n")
  gene_df = fread(parFile2, header=TRUE, data.table=FALSE) 

  stopifnot(all(colnames(valid_df) %in% colnames(gene_df)))
  
  gene_df<-gene_df[,colnames(valid_df)]

  te_df = valid_df |> dplyr::filter(Feature_gene_biotype=="transposable elements")
  
  valid_df = rbind(te_df, gene_df)
}

write.table(valid_df, file=paste0(outFile, ".count"), quote=FALSE, sep="\t", row.names=FALSE)
