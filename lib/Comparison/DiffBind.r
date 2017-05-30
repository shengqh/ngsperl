#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
pvalue<-0.05
foldChange<-1.5

library(DiffBind)

DEBUG=0
if(!DEBUG){
  configFile=args[1]
  comparisonFile=args[2]
  outputPrefix=args[3]
}else{
  configFile="/scratch/shavertm/20170512_atac-seq/bwa_macs2callpeak_narrow_CallAsSingleEnd_diffbind/result/Genotypes/Genotypes.config.txt"
  comparisonFile="/scratch/shavertm/20170512_atac-seq/bwa_macs2callpeak_narrow_CallAsSingleEnd_diffbind/result/Genotypes/Genotypes.comparison.txt"
  outputPrefix="/scratch/shavertm/20170512_atac-seq/bwa_macs2callpeak_narrow_CallAsSingleEnd_diffbind/result/Genotypes/Genotypes"
}

cat("configFile=", configFile, "\n")
cat("comparisonFile=", comparisonFile, "\n")
cat("outputPrefix=", outputPrefix, "\n")

samplesheet<-read.table(configFile, sep="\t", header=T, stringsAsFactors=F)
mb1<-dba(sampleSheet=samplesheet, bCorPlot=F)
mb1<-dba.count(mb1,score=DBA_SCORE_READS)

comparisons<-read.table(comparisonFile, se="\t", header=T, stringsAsFactors = F)
allres<-NULL
allres_sig<-NULL
idx<-1
for (idx in c(1:nrow(comparisons))){
  compName<-comparisons[idx, 1]
  compPrefix<-paste0(outputPrefix, ".", compName)
  group1name<-comparisons[idx, 2]
  group2name<-comparisons[idx, 3]
  group1<-samplesheet$Condition == group1name
  group2<-samplesheet$Condition == group2name
  mb2<-dba.contrast(mb1,group1=group1, group2=group2, name1=group1name, name2=group2name, categories=c(DBA_CONDITION), minMembers=2)
  mb2<-dba.analyze(mb2, bSubControl=FALSE, bFullLibrarySize=TRUE, bTagwise=FALSE, bCorPlot=FALSE)
  res<-dba.report(mb2,bCounts=TRUE,th=1)
  write.table(as.data.frame(res,row.names=NULL),file=paste0(compPrefix, ".tsv"),quote=F,sep="\t",row.names=F)
  select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & ((res$Fold >= foldChange) | (res$Fold <= foldChange))
  res_sig<-res[select,]
  write.table(as.data.frame(res_sig,row.names=NULL),file=paste0(compPrefix, ".sig.tsv"),quote=F,sep="\t",row.names=F)
}
