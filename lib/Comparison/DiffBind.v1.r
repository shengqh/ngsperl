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
  setwd('/scratch/cqs/shengq2/bugfix/diffbind_oldversion/result/ACBI1/')
  configFile="ACBI1.config.txt"
  comparisonFile="ACBI1.comparison.txt"
  overlapDef="ACBI1.minoverlap.txt"
  outputPrefix="ACBI1"
}

cat("configFile=", configFile, "\n")
cat("comparisonFile=", comparisonFile, "\n")
cat("overlapDef=", overlapDef, "\n")
cat("outputPrefix=", outputPrefix, "\n")

samplesheet<-read.table(configFile, sep="\t", header=T, stringsAsFactors=F)
mb1<-dba(sampleSheet=samplesheet)
#minOverlap<-ifelse(nrow(mb1) == 1, 1, 2)
#mb1<-dba.count(mb1,score=DBA_SCORE_READS, minOverlap=minOverlap)

#get consensus peaks and plot Venn for each condition, and generate count table for combined peaks across conditions 
mb1_consensus <- dba(sampleSheet=samplesheet)
overlapsheet <- read.table(overlapDef, sep="\t", header=T) 

cindex <- 1
for(cindex in c(1:nrow(overlapsheet))){
  condition <- overlapsheet[cindex,1]
  minOverlap = overlapsheet[cindex,2]
  cat("condition = ", condition, ", overlap = ", minOverlap, "\n")

  pdf(paste0(condition,".pdf",sep=""))
  dba.plotVenn(mb1,mb1$masks[names(mb1$masks)==condition][[1]])
  dev.off()
  
  masks = mb1$masks[names(mb1$masks)==condition][[1]]
  cat("  masks = ", masks, "\n")

  c_peak <- dba.peakset(mb1, masks, minOverlap=minOverlap, bRetrieve=TRUE)
  cat("  number of peaks = ", length(c_peak), "\n")
  dc_peak <- data.frame(seqnames=seqnames(c_peak),starts=start(c_peak),ends=end(c_peak),strands=strand(c_peak))
  write.table(dc_peak,file=paste0(condition,".peaks"),quote=F,row.names=F,col.names=F,sep="\t")

  mb1_consensus <- dba.peakset(mb1_consensus, masks, sampID=as.character(condition), minOverlap=minOverlap)
}

#mb1
#mb1_consensus

pdf("Overlap-condition.pdf")
dba.plotVenn(mb1_consensus, mb1_consensus$masks$Consensus)  
dev.off()

mb1_consensus <- dba(mb1_consensus,mask=mb1_consensus$masks$Consensus,minOverlap=1)
#mb1_consensus
consensus_peaks <- dba.peakset(mb1_consensus, bRetrieve=TRUE)
write.csv(consensus_peaks, "consensus_peaks.csv")
cat("Number of consensus peaks = ", length(consensus_peaks), "\n")

mb1 <- dba.count(mb1,score=DBA_SCORE_READS, peaks=consensus_peaks, bRemoveDuplicates=TRUE)
mbt = dba(mb1,bSummarizedExperiment=TRUE)
#write.table(cbind(data.frame(seqnames=seqnames(consensus_peaks),starts=start(consensus_peaks),ends=end(consensus_peaks),strands=strand(consensus_peaks)),assay(mbt)),file=paste0("counttable_raw-",outputPrefix,".txt"),quote=F,sep="\t",row.names=F)
write.table(data.frame(as.data.frame(rowRanges(mbt)),assay(mbt)),file=paste0("counttable_raw-",outputPrefix,".txt"),quote=F,sep="\t",row.names=F)
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

 #all normalizing parameters in dba.normalize
  mb2<-dba.analyze(mb2)
#  res<-dba.report(mb2,bCounts=TRUE,th=1)
  res<-dba.report(mb2,bCounts=TRUE,bNormalized=TRUE,th=1)
  write.table(as.data.frame(res,row.names=NULL),file=paste0(compPrefix, ".tsv"),quote=F,sep="\t",row.names=F)
#  select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & ((res$Fold >= foldChange) | (res$Fold <= foldChange))
  select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & (abs(res$Fold) >= log2(foldChange))
  res_sig<-res[select,]
  write.table(as.data.frame(res_sig,row.names=NULL),file=paste0(compPrefix, ".sig.tsv"),quote=F,sep="\t",row.names=F)
  write.table(as.data.frame(res_sig,row.names=NULL)[,c(1:3)],file=paste0(compPrefix, ".sig.bed"),quote=F,sep="\t",row.names=F, col.names=F)
}
