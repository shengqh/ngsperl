
###############################################
#parameters
###############################################
mafFileAllSamples=paste0(outFile,".filter.allSamples.maf")
allMafFiles=read.delim("fileList1.txt",header=FALSE,as.is=TRUE)
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
  mafFilter=filterMaf(mafFile)
  mafFilterAll=rbind(mafFilterAll,mafFilter)
}
write.table(mafFilterAll,mafFileAllSamples,sep="\t",quote = FALSE,row.names=FALSE)


###############################################
#Run mafreport on combined result
###############################################
reportOutDir=getwd()

dataForReport=initialize_maf_report_parameter(mafFileAllSamples,reportOutDir=reportOutDir,genome="hg38",
                                              reportModules=c("Initialize.Rmd","SummaryTables.Rmd","VariantsVisualization.Rmd"),
                                              dndscv.refdb="/scratch/cqs_share/references/dndscv/RefCDS_human_GRCh38.p12.rda"
                                              #                                               clinicalData="/gpfs23/scratch/cqs/zhaos/Pierre/WES/wes_sample_grp.txt",
                                              #                                               clinicalFeatures=c("PatientID","SampleGroup")
                                              #		reportTemplate="d:/source/mafreport/inst/templates/report.Rmd"
)
#mafTable=read.delim(mafFile,header=T,as.is=T,comment="#")
#mafTable=mafFilterAll
dataForReport[["maf"]]=read.maf(maf =mafFilterAll,clinicalData=dataForReport[["clinicalData"]],verbose=FALSE)

dataForReport=make_maf_report(dataForReport)

