library(data.table)
library(maftools)

args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
outputFile = args[2]
refBuild = args[3]

annovarFinalTable=fread(sep="\t",cmd=paste0("grep -v '^#' ", inputFile))
tempFile<-paste0(outputFile, ".tmp")
write.table(annovarFinalTable,tempFile,sep="\t",row.names=FALSE)
testMafTable=annovarToMaf(tempFile, refBuild=refBuild)
class(testMafTable[["AAChange"]])="character"
aaChange<-testMafTable[["AAChange"]]
aaChange[aaChange=="character(0)"]<-""
testMafTable[["AAChange"]]<-aaChange
write.table(testMafTable,outputFile,sep="\t",quote=FALSE,row.names=FALSE)
file.remove(tempFile)
