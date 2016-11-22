resultFile<-outFile
countInFastQcPreTrim<-parFile1
countInFastQcPostTrim<-parFile2
countInFastQcPostRemove<-parFile3
#resultFile<-commandArgs()[7]
#countInFastQcPreTrim<-commandArgs()[8]
#countInFastQcPostTrim<-commandArgs()[9]
#countInFastQcPostRemove<-commandArgs()[10]

library(ggplot2)
library(reshape2)
library(cowplot)

#countInFastQcPreTrim<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_pre_trim/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#countInFastQcPostTrim<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_post_trim/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#countInFastQcPostRemove<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_post_remove/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#resultFile<-"test"

if (countInFastQcPostRemove=="") {
	count1<-read.delim(countInFastQcPreTrim,header=T)
	count2<-read.delim(countInFastQcPostTrim,header=T)
	
	count1$Reads<-count1$Reads-count2$Reads
	
	count1$Label="Removed by Trimming or Removing Sequence"
	count2$Label="Reads for Mapping"
	
	countForFigure<-rbind(count2,count1)
	pdf(paste0(resultFile,".pdf"),width=14)
	print(ggplot(countForFigure,aes(x=Sample,y=Reads,fill=Label))+geom_bar(stat="identity", width=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
	dev.off()
} else {
	count1<-read.delim(countInFastQcPreTrim,header=T)
	count2<-read.delim(countInFastQcPostTrim,header=T)
	count3<-read.delim(countInFastQcPostRemove,header=T)
	
	count1$Reads<-count1$Reads-count2$Reads
	count2$Reads<-count2$Reads-count3$Reads
	
	count1$Label="Removed by Removing Sequence"
	count2$Label="Removed by Trimming"
	count3$Label="Reads for Mapping"
	
	countForFigure<-rbind(count3,count1,count2)
	pdf(paste0(resultFile,".pdf"),width=14)
	print(ggplot(countForFigure,aes(x=Sample,y=Reads,fill=Label))+geom_bar(stat="identity", width=0.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1)))
	dev.off()
}

countForFigureOut<-acast(countForFigure,Label~Sample,value.var="Reads")
write.csv(countForFigureOut,paste0(resultFile,".Reads.csv"))

