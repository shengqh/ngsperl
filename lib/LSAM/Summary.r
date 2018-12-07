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

getci=function(x)
{
  x=as.numeric(x)[2:11]
  paste0(round(mean(x),1)," (",round(mean(x)-1.96*sd(x),1),", ",round(mean(x)+1.96*sd(x),1),")")
}

T2<-NULL
for (idx in c(1:length(out1))){
  T2.idx=data.frame(LocId=out1[[idx]]$Tx.Total.DSA$Table.Tx.DSA$PatLocId,CI=apply(out1[[idx]]$Tx.Total.DSA$Table.Tx.DSA,1,getci))
  colnames(T2.idx)<-c("LocId", attr(out1[[idx]], "simName"))
  if(is.null(T2)){
    T2<-T2.idx
  }else{
    T2<-merge(T2,T2.idx,by="LocId")
  }
}

write.csv(T2,file=paste0(tablename, "_transplantTo.csv"))
