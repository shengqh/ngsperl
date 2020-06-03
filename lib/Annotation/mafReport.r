library(mafreport)

mafFileList = parSampleFile1
#reportOutDir = "."
reportOutDir = getwd()

if(parFile1 != ''){
  clinicalData = parFile1
  #need also define clinicalFeatures in rCode
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

