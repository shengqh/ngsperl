library(data.table)
library(edgeR)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    inputFile<-"/nobackup/h_cqs/nadia_sutton_projects/20231122_sutton_rnaseq_mm10/genetable/result/P2530_3428.proteincoding.count"
    outputPrefix<-'/nobackup/h_cqs/nadia_sutton_projects/20231122_sutton_rnaseq_mm10/genetable/result/P2530_3428.proteincoding'
    remove_chrM_genes<-0
  }else{
    inputFile<-args[1]
    outputPrefix<-args[2]
    if(length(args) > 2){
      remove_chrM_genes<-as.numeric(args[3])
    }else{
      remove_chrM_genes<-0
    }
  }
}

if(remove_chrM_genes){
  cat("Remove chrM genes...\n")
  with_chrM_File<-gsub(".count$", ".with_chrM.count", inputFile)
  #copy inputFile to with_chrM_File
  system(paste("cp", inputFile, with_chrM_File))

  counts<-fread(file=with_chrM_File, data.table=FALSE, check.names=FALSE)
  counts<-counts[counts$Feature_chr != "chrM", ]
  write.table(counts, file=inputFile, quote=FALSE, sep="\t", row.names=FALSE)
}

cat(inputFile, "\n")
counts<-data.frame(fread(file=inputFile), row.names=1)

name_index<-which("Feature_gene_name"==colnames(counts))
sample_index<-name_index+1
nsample=ncol(counts)
rawcounts<-counts[, c(sample_index:nsample)]
gnames=counts[, name_index, drop=FALSE]
colnames(gnames)="Gene"

#transcript level
dge=DGEList(counts=rawcounts)
dge=calcNormFactors(dge)
cpm_df=cpm(dge, normalized.lib.sizes = TRUE, log=FALSE)
stopifnot(all(rownames(gnames) == rownames(cpm_df)))

cpm_df=cbind(gnames, cpm_df)
write.csv(cpm_df, file=paste0(outputPrefix, ".cpm.csv"))

logcpm_df=cpm(dge, normalized.lib.sizes = TRUE, log=TRUE, prior.count = 2)
write.csv(logcpm_df, file=paste0(outputPrefix, ".logcpm.csv"))

#gene level
symbol_counts<-data.frame(aggregate(rawcounts, by=list(Gene=counts$Feature_gene_name), FUN=sum), row.names=1)

dge=DGEList(counts=symbol_counts)
dge=calcNormFactors(dge)
cpm_df=cpm(dge, normalized.lib.sizes = TRUE, log=FALSE)
write.csv(cpm_df, file=paste0(outputPrefix, ".gene_symbol.cpm.csv"))

logcpm_df=cpm(dge, normalized.lib.sizes = TRUE, log=TRUE, prior.count = 2)
write.csv(logcpm_df, file=paste0(outputPrefix, ".gene_symbol.logcpm.csv"))
