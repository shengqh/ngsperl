library(data.table)
library(edgeR)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    inputFile<-"/nobackup/shah_lab/shengq2/20250520_Serpina_rnaseq_hg38/genetable/result/Serpina.count"
    outputPrefix<-'/nobackup/shah_lab/shengq2/20250520_Serpina_rnaseq_hg38/genetable/result/Serpina'
    remove_chrM_genes<-1
    round_counts<-1
  }else{
    inputFile<-args[1]
    outputPrefix<-args[2]
    remove_chrM_genes<-args[3]=='1'
    round_counts<-args[4]=='1' 
  }
}

if(remove_chrM_genes | round_counts){
  suffix = ""
  if(remove_chrM_genes){
    suffix = paste0(suffix, ".remove_chrM")
  }
  if(round_counts){
    suffix = paste0(suffix, ".raw")
  }
  suffix = paste0(suffix, ".count")
  raw_File<-gsub(".count$", suffix, inputFile)

  cat("Backup", raw_File, "...\n")
  system(paste("cp", inputFile, raw_File))

  counts<-fread(file=raw_File, data.table=FALSE, check.names=FALSE)
  if(remove_chrM_genes){
    counts<-counts[counts$Feature_chr != "chrM", ]
  }
  if(round_counts){
    name_index<-which("Feature_gene_name"==colnames(counts))
    sample_index<-name_index+1
    nsample=ncol(counts)

    counts[, sample_index:nsample]<-round(counts[, sample_index:nsample])
    zero_counts<-rowSums(counts[, sample_index:nsample])==0
    if(sum(zero_counts)>0){
      cat("Remove", sum(zero_counts), "genes with all zero counts after rounding.\n")
      counts<-counts[!zero_counts, ]
    }
  }
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
