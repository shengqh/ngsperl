options(bitmapType='cairo')
options(expressions=102400)

args = commandArgs(trailingOnly = TRUE)

library(ChIPQC)
library(BiocParallel)

if(length(args) > 0){
  configFile=args[1]
  annotationName=args[2]
  chromosomes=args[3]
}else{
  configFile=r"(C:\projects\jonathan_brown\20210321_cutrun_6048_human\macs2callpeak_narrow_chipqc\result\cutrun_6048.config.txt)"
  annotationName="hg38"
  chromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
}

if (annotationName == "hg38") {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}

cat("configFile=", configFile, "\n")
cat("annotationName=", annotationName, "\n")
cat("chromosomes=", chromosomes, "\n")

register(SerialParam())

rdatafile = paste0(configFile, ".rdata")
#if (file.exists(rdatafile)){
#  load(rdatafile)
#}else{
  if(!is.na(chromosomes)){
    chromosomes = unlist(strsplit(chromosomes, split=','))
  }else{
    chromosomes = NULL
  } 

  experiment <- read.table(configFile, sep="\t", header=T)

  if(annotationName == "unknown"){
    qcresult = ChIPQC(experiment, consensus=TRUE, chromosomes=chromosomes)
  }else{
    qcresult = ChIPQC(experiment, consensus=TRUE, annotation = annotationName, chromosomes=chromosomes)
  }
  save(qcresult, file=rdatafile)
#}

ChIPQCreport(qcresult)

if(file.exists("ChIPQCreport/GenomicFeatureEnrichment.png")){
  samples=qcresult@Samples
  height=min(2000, length(names(samples)) * 500 + 1000)
  png("ChIPQCreport/GenomicFeatureEnrichment.png", width=3000, height=height, res=300)
  g<-plotRegi(qcresult)
  print(g)
  dev.off()
}
