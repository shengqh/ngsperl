options(bitmapType='cairo')
options(expressions=102400)

library(ChIPQC)
library(BiocParallel)

options_table = read.table("fileList1.txt", sep="\t")
myoptions = split(options_table$V1, options_table$V2)
annotationName=myoptions$genome
consensus=myoptions$consensus == "1"
chromosomes=myoptions$chromosomes
pcaAttributes=myoptions$pcaAttributes
pcaLabels=myoptions$pcaLabels

config_table = read.table("fileList2.txt", sep="\t")
configFile = config_table$V1[1]

if (annotationName == "hg38") {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}

cat("configFile=", configFile, "\n")
cat("annotationName=", annotationName, "\n")
cat("chromosomes=", chromosomes, "\n")
cat("consensus=", consensus, "\n")
cat("pcaAttributes=", pcaAttributes, "\n")
cat("pcaLabels=", pcaLabels, "\n")

register(SerialParam())

rdatafile = paste0(configFile, ".rdata")
if(!is.na(chromosomes)){
  chromosomes = unlist(strsplit(chromosomes, split=','))
}else{
  chromosomes = NULL
} 

experiment <- read.table(configFile, sep="\t", header=T)

if(annotationName == "unknown"){
  qcresult = ChIPQC(experiment, consensus=consensus, chromosomes=chromosomes)
}else{
  qcresult = ChIPQC(experiment, consensus=consensus, chromosomes=chromosomes, annotation=annotationName)
}
save(qcresult, file=rdatafile)

#load(rdatafile)

chipqc_version<-paste0("ChIPQC,v", packageVersion("ChIPQC"))
writeLines(chipqc_version, paste0(configFile,".ChIPQC.version"))

ChIPQCreport(qcresult)

check_unique<-function(pcaAttributes, qcresult){
  result = c()

  pcaAttributes = unlist(strsplit(pcaAttributes, ','))
  for(attr in pcaAttributes){
    if(length(unique(qcresult@DBA$samples[,attr])) > 1){
      result = c(result, attr)
    }
  }
  if(is.null(result)){
    result="ID"
  }
  return(result)
}

if(pcaAttributes != "Tissue,Factor" | pcaLabels != "Replicate"){
  pcaAttributes = check_unique(pcaAttributes, qcresult)
  pcaLabels = check_unique(pcaLabels, qcresult)
  png(file.path("ChIPQCreport","PeakPCA.png"),width=3000,height=2500,res=300)
  plotPrincomp(qcresult,attributes=pcaAttributes,label=pcaLabels)
  dev.off()
}
