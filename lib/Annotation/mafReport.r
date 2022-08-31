rm(list=ls()) 
outFile=''
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/home/shengq2/program/projects/breast_cancer_spore/20220713_wgs/BRE15136_Summary_21Oct2021_correctedv5.csv'
parFile2=''
parFile3=''
clinicalFeatures='ARM,CBR_6,TNBCtype_4,TIL_call,PDL1_IHC_Agg,BMI_CLASS,diabetes';genome='hg38'

setwd('/scratch/cqs/breast_cancer_spore/analysis/all/gatk4_13_report/result')

### Parameter setting end ###

library(mafreport)

#https://github.com/PoisonAlien/maftools/issues/532
library("wordcloud")
devtools::source_gist(id = "https://gist.github.com/PoisonAlien/3f8752a89c1d63f64afe55b441621223")

mafFileList = parSampleFile1
#reportOutDir = "."
reportOutDir = getwd()

if(parFile1 != ''){
  #need also define clinicalFeatures in rCode
  clinicalData = parFile1
  clinicalFeatures = unlist(strsplit(clinicalFeatures, ","))
  cdata<-read.csv(clinicalData, header=T, stringsAsFactor=F)
  missedFeatures<-clinicalFeatures[!(clinicalFeatures %in% colnames(cdata))]
  if (length(missedFeatures) > 0){
    stop(paste0("missing clinical features: ", paste0(missedFeatures, collapse = ",")))
  }
}else{
  clinicalData = NULL
  clinicalFeatures = NULL
}


if (!exists("genome")) {
  genome="hg19"
}

interestedGenes = NULL
if(exists("interestedGeneStr")) {
  if((!is.null(interestedGeneStr)) & (interestedGeneStr!="")) {
    interestedGeneStr = gsub("\\s+", ",", interestedGeneStr)
    interestedGenes = unlist(strsplit(interestedGeneStr, ","))
  }
}

mafFiles = read.table(mafFileList, sep="\t", header=F, stringsAsFactor=F)
for(mafFile in mafFiles$V1){
  dataForReport1=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
    vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                "Nonsense_Mutation",      "Nonstop_Mutation", "Missense_Mutation"),
    genome=genome,
    interestedGenes=interestedGenes,
    clinicalData=clinicalData,
    clinicalFeatures=clinicalFeatures
  )

  dataForReport1=make_maf_report(dataForReport1)
}

