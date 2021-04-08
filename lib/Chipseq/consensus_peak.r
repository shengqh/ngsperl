args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0){
  outFile="/scratch/cqs/shengq2/temp/concensus_peak/concensus_peak.bed"
  inputFiles=c("/scratch/jbrown_lab/shengq2/projects/20210329_cutrun_6048_human/macs2callpeak_narrow/result/BRD4/BRD4_peaks.narrowPeak.bed.gz",
               "/scratch/jbrown_lab/shengq2/projects/20210329_cutrun_6048_human/macs2callpeak_narrow/result/IgG/IgG_peaks.narrowPeak.bed")
}else{
  outFile=args[1]
  inputFiles=args[2:length(args)]
}

if(length(inputFiles) < 1 | length(inputFiles) > 5){
  stop("The number of input bed files should be between 1 to 5.")
}

options(bitmapType='cairo')

library(ChIPpeakAnno)
library(rtracklayer)

grlist=lapply(inputFiles, function(inputFile){
  toGRanges(inputFile, format="BED", header=FALSE) 
})

cat("Finding overlap peaks ...")
if(length(grlist) == 1){
  op=grlist[[1]]
} else {
  if(length(grlist) == 2){
    ol=findOverlapsOfPeaks(grlist[[1]], grlist[[2]])
  }else if(length(grlist) == 3){
    ol=findOverlapsOfPeaks(grlist[[1]], grlist[[2]], grlist[[3]])
  }else if(length(grlist) == 4){
    ol=findOverlapsOfPeaks(grlist[[1]], grlist[[2]], grlist[[3]], grlist[[4]])
  }else if(length(grlist) == 5){
    ol=findOverlapsOfPeaks(grlist[[1]], grlist[[2]], grlist[[3]], grlist[[4]], grlist[[5]])
  }
  op = ol$mergedPeaks
  op@ranges@NAMES<-paste0("ConsensusPeak_", c(1:length(op@ranges@width)))
}

export.bed(op,con=outFile)

