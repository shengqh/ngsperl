
library(reshape2)

DEBUG=1

args = commandArgs(trailingOnly=TRUE)
rnaseqFile = args[1]
snpFamFile = args[2]
outputFile = args[3]

if(DEBUG){
  rnaseqFile = "/scratch/cqs/shengq2/guoyan/prepareRnaseq/result/BRCA/BRCA.rnaseq2.count.tsv"
  snpFamFile = "/workspace/guoy1/TCGA_affy/Affy60_BC_HG19/Affy60_BC_remove_duplicate_samples.fam"
  outputFile = "/scratch/cqs/shengq2/guoyan/BRCA.samples.tsv"
}

getMeta<-function(sampleNames,analysisType){
  patients<-gsub("\\-\\d+$", "", sampleNames)
  sampleTypes<-gsub("^.+\\-","", sampleNames)
  sampleTypes[sampleTypes=="01"]<-paste0(analysisType, " Primary Solid Tumor")
  sampleTypes[sampleTypes=="06"]<-paste0(analysisType, " Metastatic")
  sampleTypes[sampleTypes=="10"]<-paste0(analysisType, " Blood Derived Normal")
  sampleTypes[sampleTypes=="11"]<-paste0(analysisType, " Solid Tissue Normal")
  meta<-data.frame(Patient=patients, SampleType=sampleTypes)
  return(meta)
}

rnaseq<-read.table(rnaseqFile, row.names=1, header=T, nrows=1, sep="\t", check.names=F)
rnaseqSamples<-colnames(rnaseq)[grepl("TCGA", colnames(rnaseq))]
rnaseqMeta<-getMeta(rnaseqSamples, "RNASeq")

snp<-read.table(snpFamFile, header=F, sep=" ", check.names=F)
snpSamples<-snp$V1
snpMeta<-getMeta(snpSamples, "SNP")

allmeta<-rbind(rnaseqMeta, snpMeta)
dmeta<-dcast(allmeta, Patient~SampleType)

write.table(dmeta, file=outputFile, sep="\t", na="")
