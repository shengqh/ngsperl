options(bitmapType='cairo')

resultFile<-outFile
countInFastQcPreTrim<-parFile1
countInFastQcTrim1<-parFile2
countInFastQcTrim2<-parFile3
#resultFile<-commandArgs()[7]
#countInFastQcPreTrim<-commandArgs()[8]
#countInFastQcTrim1<-commandArgs()[9]
#countInFastQcTrim2<-commandArgs()[10]

library(ggplot2)
library(reshape2)
library(cowplot)

#countInFastQcPreTrim<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_pre_trim/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#countInFastQcTrim1<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_post_trim/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#countInFastQcTrim2<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_post_remove/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#resultFile<-"test"
readCount<-function(fileName){
  count1<-read.delim(fileName,header=T,check.names=F)
  count1<-count1[!duplicated(count1$Sample), ]
  rownames(count1)<-count1$Sample
  return(count1)
}

count1<-readCount(countInFastQcPreTrim)
if (countInFastQcTrim1=="") {
  count1$Label="Reads for Mapping"
  countForFigure<-count1
} else {
if (countInFastQcTrim2=="") {
  count2<-readCount(countInFastQcTrim1)
  
  count1<-count1[rownames(count2),]
  count1$Reads<-count1$Reads-count2$Reads
  
  if(grepl("post_remove", countInFastQcTrim1)){
    count1$Label="Removed by Removing Sequence"
  } else if (grepl("post_trim", countInFastQcTrim1)){
    count1$Label="Removed by Trimming"
  } else {
    count1$Label="Removed by Trimming or Removing Sequence"
  }  
  
  count2$Label="Reads for Mapping"
  
  countForFigure<-rbind(count2,count1)
} else {
  count2<-readCount(countInFastQcTrim1)
  count3<-readCount(countInFastQcTrim2)
  
  count1<-count1[rownames(count3),]
  count2<-count2[rownames(count3),]
  
  count1$Reads<-count1$Reads-count2$Reads
  count2$Reads<-count2$Reads-count3$Reads
  
  count1$Label="Removed by Removing Sequence"
  count2$Label="Removed by Trimming"
  count3$Label="Reads for Mapping"
  
  countForFigure<-rbind(count3,count1,count2)
}
}

width<-max(14, length(unique(countForFigure$Sample)) * 0.15)
pdf(paste0(resultFile,".pdf"),width=width)
print(ggplot(countForFigure,aes(x=Sample,y=Reads,fill=Label))+geom_bar(stat="identity", width=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
dev.off()

countForFigureOut<-acast(countForFigure,Label~Sample,value.var="Reads")
write.csv(countForFigureOut,paste0(resultFile,".Reads.csv"))

