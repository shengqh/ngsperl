options(bitmapType='cairo')

args = commandArgs(trailingOnly = TRUE)
pvalue<-0.05
foldChange<-1.5

library(DiffBind)

DEBUG=0
if(!DEBUG){
  configFile=args[1]
  comparisonFile=args[2]

# the new parameter defining minoverlap for each condition 
  overlapDef=args[3]
  outputPrefix=args[4]
}else{
  configFile="/scratch/shavertm/20170512_atac-seq/bwa_macs2callpeak_narrow_CallAsSingleEnd_diffbind/result/Genotypes/Genotypes.config.txt"
  comparisonFile="/scratch/shavertm/20170512_atac-seq/bwa_macs2callpeak_narrow_CallAsSingleEnd_diffbind/result/Genotypes/Genotypes.comparison.txt"
  outputPrefix="/scratch/shavertm/20170512_atac-seq/bwa_macs2callpeak_narrow_CallAsSingleEnd_diffbind/result/Genotypes/Genotypes"
}

cat("configFile=", configFile, "\n")
cat("comparisonFile=", comparisonFile, "\n")
cat("overlapDef=", overlapDef, "\n")
cat("outputPrefix=", outputPrefix, "\n")

samplesheet<-read.table(configFile, sep="\t", header=T, stringsAsFactors=F)
mb1<-dba(sampleSheet=samplesheet, bCorPlot=F)
#minOverlap<-ifelse(nrow(mb1) == 1, 1, 2)
#mb1<-dba.count(mb1,score=DBA_SCORE_READS, minOverlap=minOverlap)

#get consensus peaks and plot Venn for each condition, and generate count table for combined peaks across conditions 
mb1_consensus <- dba(sampleSheet=samplesheet, bCorPlot=F)
overlapsheet <- read.table(overlapDef, sep="\t", header=T) 

cindex <- 1
for(cindex in c(1:nrow(overlapsheet))){
  condition <- overlapsheet[cindex,1]
  pdf(paste0(condition,".pdf",sep=""))
  dba.plotVenn(mb1,mb1$masks[names(mb1$masks)==condition][[1]])
  dev.off()
  
  c_peak <- dba.peakset(mb1, mb1$masks[names(mb1$masks)==condition][[1]], minOverlap=overlapsheet[cindex,2], bRetrieve=TRUE)
  dc_peak <- data.frame(seqnames=seqnames(c_peak),starts=start(c_peak),ends=end(c_peak),strands=strand(c_peak))
  write.table(dc_peak,file=paste0(condition,".peaks"),quote=F,row.names=F,col.names=F,sep="\t")
  mb1_consensus <- dba.peakset(mb1_consensus, mb1_consensus$masks[names(mb1_consensus$masks)==condition][[1]], sampID=as.character(condition), minOverlap=overlapsheet[cindex,2])
}

#mb1
#mb1_consensus

pdf("Overlap-condition.pdf")
dba.plotVenn(mb1_consensus, mb1_consensus$masks$Consensus)  
dev.off()

mb1_consensus <- dba(mb1_consensus,mask=mb1_consensus$masks$Consensus,minOverlap=1)
#mb1_consensus
consensus_peaks <- dba.peakset(mb1_consensus, bRetrieve=TRUE)
mb1 <- dba.count(mb1,score=DBA_SCORE_READS, peaks=consensus_peaks, bRemoveDuplicates=TRUE)
mbt = dba(mb1,bSummarizedExperiment=TRUE)
write.table(cbind(data.frame(seqnames=seqnames(consensus_peaks),starts=start(consensus_peaks),ends=end(consensus_peaks),strands=strand(consensus_peaks)),assay(mbt)),file=paste0("counttable_raw-",outputPrefix,".txt"),quote=F,sep="\t",row.names=F)
mb1


comparisons<-read.table(comparisonFile, sep="\t", header=T, stringsAsFactors = F)
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
#  mb2<-dba.contrast(mb1,group1=group1, group2=group2, name1=group1name, name2=group2name, categories=c(DBA_CONDITION), minMembers=2)
  mb2<-dba.contrast(mb1,group1=group1, group2=group2, name1=group1name, name2=group2name)

  mb2<-dba.analyze(mb2, bSubControl=FALSE, bFullLibrarySize=TRUE, bTagwise=FALSE, bCorPlot=FALSE)
#  res<-dba.report(mb2,bCounts=TRUE,th=1)
  res<-dba.report(mb2,bCounts=TRUE,bNormalized=TRUE,th=1)
  write.table(as.data.frame(res,row.names=NULL),file=paste0(compPrefix, ".tsv"),quote=F,sep="\t",row.names=F)
#  select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & ((res$Fold >= foldChange) | (res$Fold <= foldChange))
  select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & (abs(res$Fold) >= log2(foldChange))
  res_sig<-res[select,]
  write.table(as.data.frame(res_sig,row.names=NULL),file=paste0(compPrefix, ".sig.tsv"),quote=F,sep="\t",row.names=F)
}
