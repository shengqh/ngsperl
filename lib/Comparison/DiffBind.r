#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
pvalue<-0.05
foldChange<-1.5

library(DiffBind)

DEBUG=0
if(!DEBUG){
  configFile=args[1]
  outputPrefix=args[2]
}else{
  configFile="/scratch/cqs/shengq1/brown/20170109_atacseq_3585_13to19/bwa_macs2callpeak_broad_diffbind/result/Dex_vs_NoTx/Dex_vs_NoTx.txt"
  outputPrefix="/scratch/cqs/shengq1/brown/20170109_atacseq_3585_13to19/bwa_macs2callpeak_broad_diffbind/result/Dex_vs_NoTx/Dex_vs_NoTx.result"
}

cat("configFile=", configFile, "\n")
cat("outputPrefix=", outputPrefix, "\n")

samplesheet<-read.table(configFile, sep="\t", header=T)
mb1<-dba(sampleSheet=samplesheet)
mb1<-dba.count(mb1,score=DBA_SCORE_READS)
mb1<-dba.contrast(mb1,categories=DBA_CONDITION,minMembers=2)
mb1<-dba.analyze(mb1, bSubControl=FALSE, bFullLibrarySize=TRUE, bTagwise=FALSE)
res<-dba.report(mb1,bCounts=TRUE,th=1)
write.table(as.data.frame(res,row.names=NULL),file=paste0(outputPrefix,".tsv"),quote=F,sep="\t",row.names=F)
select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & ((res$Fold >= foldChange) | (res$Fold <= foldChange))
res_sig<-res[select,]
write.table(as.data.frame(res_sig,row.names=NULL),file=paste0(outputPrefix,".sig.tsv"),quote=F,sep="\t",row.names=F)
