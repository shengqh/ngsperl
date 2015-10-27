resultFile<-commandArgs()[7]
countInFastQcPreTrim<-commandArgs()[8]
countInFastQcPostTrim<-commandArgs()[9]
countInFastQcPostRemove<-commandArgs()[9]

library(ggplot2)
library(reshape)

#countInFastQcPreTrim<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_pre_trim/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#countInFastQcPostTrim<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_post_trim/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#countInFastQcPostRemove<-"/scratch/cqs/zhaos/vickers/20151017_3018-KCV-45-46/fastqc_post_remove/result/3018-KCV-45-46.FastQC.summary.reads.tsv"
#resultFile<-"test"

count1<-read.delim(countInFastQcPreTrim,header=T)
count2<-read.delim(countInFastQcPostTrim,header=T)
count3<-read.delim(countInFastQcPostRemove,header=T)

count1[,2]<-count1[,2]-count2[,2]
count2[,2]<-count2[,2]-count3[,2]

count1$Label="Removed by Trim"
count2$Label="Removed by Removing Sequence"
count3$Label="Reads for Mapping"

countForFigure<-rbind(count1,count2,count3)
pdf(paste0(resultFile,".pdf"),width=14)
ggplot(countForFigure,aes(x=Sample,y=Reads,fill=Label))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


