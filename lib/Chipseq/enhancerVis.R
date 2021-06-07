
options(bitmapType='cairo')

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
    signalData<-read.delim(signalFile, stringsAsFactor=F,check.names=F)
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
  mfold<-NULL
  mann<-NULL
  mpoint<-NULL
  for(row in c(2:(ncol(allData)-1))){
    rName<-colnames(allData)[row]
    for(col in c((row+1):ncol(allData))){
      cName<-colnames(allData)[col]
      sdata<-allData[,c(row,col)]
      sdata<-sdata[rowSums(sdata)>0,]
      
      mdata<-melt(sdata)
      mdata$variable<-as.character(mdata$variable)
      mdata$Row<-rName
      mdata$Col<-cName
      mdata$Name<-mdata$variable
      mdata$Name[mdata$variable==colnames(allData)[row]]<-"Row"
      mdata$Name[mdata$variable==colnames(allData)[col]]<-"Col"
      mtss<-rbind(mtss, mdata)
      
      colnames(sdata)<-c("Y", "X")
      sdata$Col<-rName
      sdata$Row<-cName
      mpoint<-rbind(mpoint, sdata)
      
      fdata=data.frame(Log2Fold=log2(sdata[,2] / sdata[,1]))
      fdata$Row<-rName
      fdata$Col<-cName
      mfold<-rbind(mfold, fdata)
      
      pvalue<-cor.test(sdata[,1], sdata[,2], alternative="two.sided", method="spearman")$p.value
      ann_text<-data.frame(Row=rName, Col=cName, pValue=paste0("pValue=",pvalue))
      mann<-rbind(mann, ann_text)
    }
  }
  colnames(mtss)<-c("Sample", "Signal", "Row", "Col", "Name")
  mtss$Signal<-log(mtss$Signal+1)
  
  png(file=paste0(resultFilePrefix, ".", colName, ".pair.violin.png"), width=3000, height=3000, res=300)
  g<-ggplot(mtss) + geom_violin(aes(x=Name, y=Signal, col=Sample)) + 
    facet_grid(Row~Col) +
    theme(axis.text.x = element_blank(), legend.position="top") + xlab("")
  print(g)
  dev.off()
  
  png(file=paste0(resultFilePrefix, ".", colName, ".pair.corr.png"), width=3000, height=3000, res=300)
  g<-ggplot(mpoint, aes(x=X, y=Y)) +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red")+
    facet_grid(Row~Col, switch = "both") +
    xlab("") + ylab("")
  print(g)
  dev.off()
  
  png(file=paste0(resultFilePrefix, ".", colName, ".pair.fold.png"), width=3000, height=3000, res=300)
  g<-ggplot(mfold, aes(x=1,y=Log2Fold)) + 
    geom_violin() + 
    facet_grid(Row~Col)+
    theme(axis.text.x = element_blank(), legend.position="top") + xlab("") + background_grid(major = "xy", minor = "none") + 
    ylab("log2(Column/Row)") + geom_hline(aes(yintercept=0), color="red")
  print(g)
  dev.off()
  
  write.table(mann, file=paste0(resultFilePrefix, ".", colName, ".spearman.tsv"), row.names=F, sep="\t", quote=F)
}

alldata<-NULL
i<-1
for (i in c(1:nrow(signalFileList))){
  signalFile<-signalFileList[i,1]
  signalName<-signalFileList[i,2]
  signalData<-read.delim(signalFile, stringsAsFactor=F,check.names=F)
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
