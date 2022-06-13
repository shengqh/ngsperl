
options(bitmapType='cairo')

library(reshape2)
library(ggplot2)
library(cowplot)

resultFile<-outFile
annovarFiles<-parSampleFile1

annovarFileList<-read.delim(annovarFiles,header=F,check.names=F,stringsAsFactor=F)

alldata<-NULL
i<-1
for (i in c(1:nrow(annovarFileList))){
  annovarFile<-annovarFileList[i,1]
  annovarName<-annovarFileList[i,2]
  annovarData<-read.delim(annovarFile, stringsAsFactor=F,check.names=F)
  annovarData$Sample=annovarName
  annovarData<-annovarData[,c("Chr","Func.refGene","Sample")]
  alldata<-rbind(alldata, annovarData)
}

saveRDS(alldata, "alldata.rds")
