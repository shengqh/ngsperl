library(ggplot2)
library(data.table)
library(dplyr)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    inputFile<-"/scratch/cqs/shengq2/justin_balko_projects/20220209_rnaseq_7312_mm10/genetable/result/P7312.count"
    outputPrefix<-'/scratch/cqs/shengq2/justin_balko_projects/20220209_rnaseq_7312_mm10/genetable/result/P7312.tpm.csv'
  }else{
    inputFile<-args[1]
    outputPrefix<-args[2]
  }
}

cat(inputFile, "\n")
counts<-read.table(file=inputFile, sep="\t", header=T, row.names=1)

name_index<-which("Feature_gene_name"==colnames(counts))
sample_index<-name_index+1
nsample=ncol(counts)
rawcounts<-counts[,c(sample_index:nsample)]
featureLength=counts$Feature_length
id_name_map=split(counts$Feature_gene_name, rownames(counts))
name_id_map=split(rownames(counts),counts$Feature_gene_name)

counts_to_tpm <- function(counts, featureLength) {
  stopifnot(length(featureLength) == nrow(counts))
  
  logFeatureLength<-log(featureLength)
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - logFeatureLength[i]
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

gene_to_symbol<-function(tpm, id_name_map, name_id_map){
  dup_symbols=name_id_map[lapply(name_id_map, length)>1]

  unique_ids=unlist(name_id_map[!(names(name_id_map) %in% names(dup_symbols))])
  unique_tpm=tpm[unique_ids,]
  rownames(unique_tpm)<-id_name_map[rownames(unique_tpm)]

  rm_dup_tpm=NULL
  symbol=names(dup_symbols)[1]
  for(symbol in names(dup_symbols)){
    dup_tpm=tpm[rownames(tpm) %in% unlist(dup_symbols[symbol]),]
    merged_tpm=colSums(dup_tpm)
    rm_dup_tpm<-rbind(rm_dup_tpm, merged_tpm)
  }

  rownames(rm_dup_tpm)<-names(dup_symbols)

  result<-rbind(unique_tpm, rm_dup_tpm)
  return(result)
}

tpm<-counts_to_tpm(rawcounts, featureLength)
symbol_tpm<-gene_to_symbol(tpm, id_name_map, name_id_map)
write.csv(symbol_tpm, file=outputPrefix)
