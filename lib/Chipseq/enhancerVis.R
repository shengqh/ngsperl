library(reshape2)
library(ggplot2)
library(cowplot)

resultFile<-outFile
signalFiles<-parSampleFile1

signalFileList<-read.delim(signalFiles,header=F,check.names=F,stringsAsFactor=F)

process<-function(signalFileList, colIndex, colName, resultFilePrefix){
  allData<-NULL
  i<-1
  for (i in c(1:nrow(signalFileList))){
    signalFile<-signalFileList[i,1]
    signalName<-signalFileList[i,2]
    signalData<-read.delim(signalFile, stringsAsFactor=F)
    curData<-signalData[,c(1,colIndex)]
    colnames(curData)<-c("Gene", signalName)
    
    if(is.null(allData)){
      allData<-curData
    }else{
      allData<-merge(allData, curData, by="Gene")
    }
  }
  
  allData<-allData[rowSums(allData[,-1,drop=F])>0,]
  write.table(allData, file=paste0(resultFilePrefix, ".", colName, ".tsv"), sep="\t", row.names=F, quote=F)
  
  mtss<-NULL
  for(row in c(2:(ncol(allData)-1))){
    for(col in c((row+1):ncol(allData))){
      sdata<-allData[,c(row,col)]
      sdata<-sdata[rowSums(sdata)>0,]
      mdata<-melt(sdata)
      mdata$variable<-as.character(mdata$variable)
      mdata$Row<-colnames(allData)[row]
      mdata$Col<-colnames(allData)[col]
      mdata$Name<-mdata$variable
      mdata$Name[mdata$variable==colnames(allData)[row]]<-"Row"
      mdata$Name[mdata$variable==colnames(allData)[col]]<-"Col"
      mtss<-rbind(mtss, mdata)
    }
  }
  colnames(mtss)<-c("Sample", "Signal", "Row", "Col", "Name")
  mtss$Signal<-log(mtss$Signal+1)
  
  png(file=paste0(resultFilePrefix, ".", colName, ".pair.violin.png"), width=3000, height=3000, res=300)
  g<-ggplot(mtss, aes(x=Name, y=Signal, col=Sample)) + geom_violin() + facet_grid(Row~Col)+
    theme(axis.text.x = element_blank(), legend.position="top") + xlab("")
  print(g)
  dev.off()
}

alldata<-NULL
i<-1
for (i in c(1:nrow(signalFileList))){
  signalFile<-signalFileList[i,1]
  signalName<-signalFileList[i,2]
  signalData<-read.delim(signalFile, stringsAsFactor=F)
  colnames(signalData)<-c("Gene", "TSS", "Distal")
  signalData<-signalData[rowSums(signalData[,c(-1)])>0,]
  signalMeltData<-melt(signalData, id="Gene")
  signalMeltData$Sample<-signalName
  alldata<-rbind(alldata, signalMeltData)
}

png(file=paste0(resultFile, ".boxplot.png"), width=3000, height=2000, res=300)
g<-ggplot(alldata, aes(x=Sample, y=value)) + geom_boxplot() + facet_grid(.~variable) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Signal")
print(g)
dev.off()

png(file=paste0(resultFile, ".violin.png"), width=3000, height=2000, res=300)
g<-ggplot(alldata, aes(x=Sample, y=value)) + geom_violin() + facet_grid(.~variable) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Signal")
print(g)
dev.off()

alldata$value<-log(alldata$value + 1)
png(file=paste0(resultFile, ".log.boxplot.png"), width=3000, height=2000, res=300)
g<-ggplot(alldata, aes(x=Sample, y=value)) + geom_boxplot() + facet_grid(.~variable) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Signal")
print(g)
dev.off()

png(file=paste0(resultFile, ".log.violin.png"), width=3000, height=2000, res=300)
g<-ggplot(alldata, aes(x=Sample, y=value)) + geom_violin() + facet_grid(.~variable) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Signal")
print(g)
dev.off()

process(signalFileList, 2, "tss", resultFile)
process(signalFileList, 3, "distal", resultFile)
