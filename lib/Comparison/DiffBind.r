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
  optionFile=args[5]
}else{
  setwd('/scratch/cqs/shengq2/bugfix/diffbind_minOverlap2_summits1/result/ACBI1')
  configFile="ACBI1.config.txt"
  comparisonFile="ACBI1.comparison.txt"
  overlapDef="ACBI1.minoverlap.txt"
  outputPrefix="ACBI1"
  optionFile="fileList1.txt"
}

cat("configFile=", configFile, "\n")
cat("comparisonFile=", comparisonFile, "\n")
cat("overlapDef=", overlapDef, "\n")
cat("outputPrefix=", outputPrefix, "\n")
cat("optionFile=", optionFile, "\n")

get_summit<-function(value){
  if(is.null(value)){
    return(FALSE)
  }
  if(is.na(value)){
    return(FALSE)
  }
  if(value == '0'){
    return(FALSE)
  }
  return(as.numeric(value))
}

get_consensus_minOverlap<-function(value, default_value=2){
  if(is.null(value)){
    return(default_value)
  }
  if(is.na(value)){
    return(default_value)
  }
  return(as.numeric(value))
}

option_table<-read.table(optionFile, sep="\t")
myoptions = split(option_table$V1, option_table$V2)

summits=get_summit(myoptions$summits)
cat("summits=", summits, "\n")

consensus_minOverlap = get_consensus_minOverlap(myoptions$consensus_minOverlap, 2)
cat("consensus_minOverlap=", consensus_minOverlap, "\n")

sampleSheet<-read.table(configFile, sep="\t", header=T, stringsAsFactors=F)
sampleSheet

comparisons<-read.table(comparisonFile, sep="\t", header=T, stringsAsFactors = F)
all_groups<-unique(c(comparisons$Group1, comparisons$Group2))
missed_groups<-all_groups[!all_groups %in% sampleSheet$Condition]
if(length(missed_groups) > 0){
  stop("Group defined in comparison not found in config file:", paste0(missed_groups, collapse = ","))
}

overlapsheet <- read.table(overlapDef, sep="\t", header=T) 
if(nrow(overlapsheet) > 0){
  missed_groups<-overlapsheet$Condition[!overlapsheet$Condition %in% sampleSheet$Condition]
  if(length(missed_groups) > 0){
    stop("Group defined in MinOverlap not found in config file:", paste0(missed_groups, collapse = ","))
  }
}

is_file_exists<-file.exists(sampleSheet$Peaks)
if(any(!is_file_exists)){
  not_exists=sampleSheet$Peaks[!is_file_exists]
  stop(paste0("Those files not exits: \n", paste0(not_exists, collapse = "\n")))
}

is_file_exists<-file.exists(sampleSheet$bamReads)
if(any(!is_file_exists)){
  not_exists=sampleSheet$bamReads[!is_file_exists]
  stop(paste0("Those files not exits: \n", paste0(not_exists, collapse = "\n")))
}

mb1<-dba(sampleSheet=sampleSheet, peakFormat="bed")

png(paste0(outputPrefix, ".heatmap.png"), width=2000, height=2000, res=300)
plot(mb1)
dev.off()

if(nrow(overlapsheet) > 0){
  cat("checking overlap ... \n")
  #get consensus peaks and plot Venn for each condition, and generate count table for combined peaks across conditions 
  mb1_consensus <- dba(sampleSheet=sampleSheet, minOverlap = consensus_minOverlap)

  cindex <- 1
  for(cindex in c(1:nrow(overlapsheet))){
    condition <- overlapsheet[cindex,1]
    minOverlap <- overlapsheet[cindex,2]
    cat("condition = ", condition, ", overlap = ", minOverlap, "\n")

    if (!condition %in% names(mb1$masks)){
      all_conditions = unique(sampleSheet$Condition)
      stop(paste0("condition ", condition, " in overlap file was not found in Condition of configuration file: ", paste0(all_conditions, collapse = ", ")))
    }
    
    masks = mb1$masks[names(mb1$masks)==condition][[1]]
    #cat("  masks = ", masks, "\n")
    samples = mb1$sample$SampleID[masks]
    cat("  samples = ", samples, "\n")

    png(paste0(outputPrefix, ".", condition, ".overlap.venn.png"), width=2000, height=2000, res=300)
    dba.plotVenn(mb1, masks)
    dev.off()
    
    pdf(paste0(outputPrefix, ".", condition, ".overlap.venn.pdf"))
    dba.plotVenn(mb1, masks)
    dev.off()
    
    c_peak <- dba.peakset(mb1, masks, minOverlap=minOverlap, bRetrieve=TRUE)
    cat("  number of peaks = ", length(c_peak), "\n")

    dc_peak <- data.frame(seqnames=seqnames(c_peak),starts=start(c_peak),ends=end(c_peak),strands=strand(c_peak))
    write.table(dc_peak,file=paste0(outputPrefix, ".", condition,".overlap.peaks.txt"),quote=F,row.names=F,col.names=F,sep="\t")
    
    mb1_consensus <- dba.peakset(DBA=mb1_consensus, peaks=c_peak, sampID=as.character(condition), consensus=TRUE, minOverlap=1)
  }

  png(paste0(outputPrefix, ".overlap-condition.png"), width =2000, height=2000, res=300)
  dba.plotVenn(mb1_consensus, mb1_consensus$masks$Consensus)  
  dev.off()

  pdf(paste0(outputPrefix, ".overlap-condition.pdf"))
  dba.plotVenn(mb1_consensus, mb1_consensus$masks$Consensus)  
  dev.off()

  #set minOverlap=1 between different sample group.
  mb1_consensus <- dba(mb1_consensus, mask=mb1_consensus$masks$Consensus, minOverlap=1)
  consensus_peaks <- dba.peakset(mb1_consensus, bRetrieve=TRUE)
  write.csv(consensus_peaks, paste0(outputPrefix, ".consensus_peaks.csv"))
  cat("Number of consensus peaks = ", length(consensus_peaks), "\n")

  rm(mb1_consensus)

  #now get count for all samples based on consensus_peaks.
  mb1 <- dba.count(mb1, summits=summits, score=DBA_SCORE_READS, peaks=consensus_peaks, bRemoveDuplicates=TRUE)
}else{
  #now get count for all samples
  mb1 <- dba.count(mb1, summits=summits, score=DBA_SCORE_READS, bRemoveDuplicates=TRUE)
}

mbt = dba(mb1,bSummarizedExperiment=TRUE)

write.table(data.frame(as.data.frame(rowRanges(mbt)),assay(mbt)),
            file=paste0(outputPrefix,".counttable_raw.txt"),quote=F,sep="\t",row.names=F)

allres<-NULL
allres_sig<-NULL
idx<-1
for (idx in c(1:nrow(comparisons))){
  compName<-comparisons[idx, 1]
  compPrefix<-paste0(outputPrefix, ".", compName)

  group1name<-comparisons[idx, 2]
  group2name<-comparisons[idx, 3]
  
  group1<-sampleSheet$Condition == group1name
  group2<-sampleSheet$Condition == group2name

  mb2<-dba.contrast(mb1, group1=group1, group2=group2, name1=group1name, name2=group2name)

  #all normalizing parameters in dba.normalize
  mb2<-dba.analyze(mb2)

  res<-dba.report(mb2,bCounts=TRUE,bNormalized=TRUE,th=1)
  write.table(as.data.frame(res,row.names=NULL),file=paste0(compPrefix, ".tsv"),quote=F,sep="\t",row.names=F)

  select<-(!is.na(res$FDR)) & (res$FDR<pvalue) & (abs(res$Fold) >= log2(foldChange))
  res_sig<-res[select,]

  write.table(as.data.frame(res_sig,row.names=NULL),file=paste0(compPrefix, ".sig.tsv"),quote=F,sep="\t",row.names=F)
  
  write.table(as.data.frame(res_sig,row.names=NULL)[,c(1:3)],file=paste0(compPrefix, ".sig.bed"),quote=F,sep="\t",row.names=F, col.names=F)
}
