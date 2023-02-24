options(bitmapType='cairo')
options(expressions=102400)

library(ChIPQC)
library(BiocParallel)

options_table = read.table("fileList1.txt", sep="\t")
myoptions = split(options_table$V1, options_table$V2)
annotationName=myoptions$genome
consensus=myoptions$consensus == "1"
chromosomes=myoptions$chromosomes

config_table = read.table("fileList2.txt", sep="\t")
configFile = config_table$V1[1]

if (annotationName == "hg38") {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
}

cat("configFile=", configFile, "\n")
cat("annotationName=", annotationName, "\n")
cat("chromosomes=", chromosomes, "\n")
cat("consensus=", consensus, "\n")

register(SerialParam())

rdatafile = paste0(configFile, ".rata")
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
