library(mafreport)

mafFileList = parSampleFile1
reportOutDir = "."
clinicalData = parFile1
if(!is.null(interestedGeneStr)) {
  interestedGeneStr = gsub("\\s+", ",", interestedGeneStr)
  interestedGenes = unlist(strsplit(interestedGeneStr, ","))
}else{
  interestedGenes = NULL
}

mafFiles = read.table(mafFileList, sep="\t", header=F, stringsAsFactor=F)
for(mafFile in mafFiles$V1){
  dataForReport1=initialize_maf_report_parameter(mafFile,reportOutDir=reportOutDir,
    vc_nonSyn=c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                "Nonsense_Mutation",      "Nonstop_Mutation", "Missense_Mutation"),
    interestedGenes=interestedGenes,
    clinicalData=clinicalData,
    clinicalFeatures=clinicalFeatures
  )

  dataForReport1=make_maf_report(dataForReport1)
}
