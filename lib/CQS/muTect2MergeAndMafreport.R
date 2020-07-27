
#other parameters defined by ExomeSeq.pm
#genome

###############################################
#parameters
###############################################
mafFileAllSamples=paste0(outFile,".filter.allSamples.maf")
allMafFiles=read.delim("fileList1.txt",header=FALSE,as.is=TRUE)
clinicalData=parFile1
library(mafreport)


###############################################
#filter results and merge all samples
###############################################
mafFilterAll=NULL #combined data for report
for (i in 1:nrow(allMafFiles)) {
  mafFile=allMafFiles[i,1]
  if (file.exists(mafFile)) {
    print(mafFile)
  } else {
    next;
  }
  mafFilter=filterMaf(mafFile,mafMax = 0.01,ExAC_FILTER=FALSE)
  mafFilterAll=rbind(mafFilterAll,mafFilter)
}
#write.table(mafFilterAll,mafFileAllSamples,sep="\t",quote = FALSE,row.names=FALSE)
write.table(mafFilterAll,mafFileAllSamples,sep="\t",row.names=FALSE)


###############################################
#Run mafreport on combined result
###############################################
reportOutDir=getwd()

if (genome=="hg19") {
  dndscv.refdb=genome
} else if (1) {
  dndscv.refdb="/scratch/cqs_share/references/dndscv/RefCDS_human_GRCh38.p12.rda"
}

dataForReport=initialize_maf_report_parameter(paste0(getwd(),"/",mafFileAllSamples),reportOutDir=reportOutDir,genome=genome,
                                              #reportModules=c("Initialize.Rmd","SummaryTables.Rmd","VariantsVisualization.Rmd"),
                                              dndscv.refdb=dndscv.refdb,
                                              clinicalData=clinicalData,
                                              clinicalFeatures=clinicalFeatures
                                              #		reportTemplate="d:/source/mafreport/inst/templates/report.Rmd"
)
#mafTable=read.delim(mafFile,header=T,as.is=T,comment="#")
#mafTable=mafFilterAll
#dataForReport[["maf"]]=read.maf(maf =mafFilterAll,clinicalData=dataForReport[["clinicalData"]],verbose=FALSE)

dataForReport=make_maf_report(dataForReport)

