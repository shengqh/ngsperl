options(bitmapType='cairo')

resultPrefix<-outFile
countInFastQcRaw<-parFile1
countInFastQcPostJoin<-parFile2
countMapped<-parFile3

library(ggplot2)
library(cowplot)
library(reshape2)

#resultPrefix<-"/workspace/shengq1/guoyan/20160922_rnaediting/test"
#countInFastQcRaw<-"/workspace/shengq1/guoyan/20160922_rnaediting/fastqc_raw/result/rnaediting.FastQC.summary.reads.tsv"
#countInFastQcPostJoin<-"/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_fastqc/result/rnaediting.FastQC.summary.reads.tsv"
#countMapped<-"/workspace/shengq1/guoyan/20160922_rnaediting/fastq_join_bowtie/result/rnaediting.MappedStat.Reads.tsv"

count1<-read.delim(countInFastQcRaw,header=T)
rownames(count1)=count1$Sample
count2<-read.delim(countInFastQcPostJoin,header=T)
rownames(count2)=count2$Sample
count3<-read.delim(countMapped,header=T)
rownames(count3)=count3$Sample
count3$Reads=count3$mapped
count3=count3[,c("Sample", "Reads")]

count2=count2[rownames(count1),]
count3=count3[rownames(count1),]

count1$Reads<-count1$Reads-count2$Reads
count2$Reads<-count2$Reads-count3$Reads

count1$Category="Removed"
count2$Category="Unmapped"
count3$Category="Mapped"

countForFigure<-rbind(count1,count2,count3)
countForFigure$Category=factor(countForFigure$Category, levels=c("Removed", "Unmapped", "Mapped"))

width=max(2000, 100 * nrow(count1))
png(file=paste0(resultPrefix,".Reads.png"), height=2000, width=width, res=300)
g=ggplot(countForFigure,aes(x=Sample,y=Reads,fill=Category))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        legend.position="top")
print(g)
dev.off()

countForFigureOut<-acast(countForFigure,Category~Sample,value.var="Reads")
countForFigureOut=data.frame(t(countForFigureOut))
countForFigureOut$Total=rowSums(countForFigureOut)
countForFigureOut$Sample=rownames(countForFigureOut)
countForFigureOut=countForFigureOut[,c(5,4,1,2,3)]

write.table(countForFigureOut,paste0(resultPrefix,".Reads.tsv"),row.names=F,quote=F,sep="\t")
