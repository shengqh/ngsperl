require(rms)
require(xlsx)
require(data.table)
require(Hmisc)

tablename<-outFile
ParentDataFolder<-parFile1
setwd(ParentDataFolder)

ctrFile<-parFile2

load(ctrFile)
options("stringsAsFactors"=FALSE)
out1=OutputSummary(ParentDataFolder,ctr.opo)
GenerateTables(out1,csv=T,tablename=tablename)