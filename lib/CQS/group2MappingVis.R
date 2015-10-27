resultFile<-commandArgs()[7]
bacCountFile<-commandArgs()[8]
bacNameFile<-commandArgs()[9]

library(ggplot2)
library(reshape)

#setwd("/gpfs21/scratch/cqs/zhaos/vickers/20150911_3018-KCV-17-18-nonHuman")
#bacCountFile<-"./bowtie1_bacteria_group2_pm_table/result/bacteria_group2_pm_3018-KCV-17-18-nonHuman.count"
#bacNameFile<-"/scratch/cqs/zhaos/vickers/reference/bacteria/group2/20150902.log"

mappedCount<-read.delim(bacCountFile,header=T,row.names=1,as.is=T)
bacName<-read.delim(bacNameFile,header=F,skip=1,row.names=2,as.is=T)

row.names(mappedCount)<-gsub("\\.\\d$","",row.names(mappedCount))
row.names(bacName)<-gsub("\\.fna$","",row.names(bacName))
row.names(mappedCount)<-bacName[row.names(mappedCount),1]
write.csv(mappedCount,paste0(bacCountFile,".NameToGenome.csv"))

#barplot(mappedCount)
mappedCountForFigure<-mappedCount
mappedCountForFigure$Bacteria<-row.names(mappedCountForFigure)
mappedCountForFigure<-melt(mappedCountForFigure)
colnames(mappedCountForFigure)<-c("Bacteria","Sample","Reads")
mappedCountForFigure$Sample<-gsub("^X","",mappedCountForFigure$Sample)

pdf(paste0(resultFile,".pdf"),width=14)
ggplot(mappedCountForFigure,aes(x=Sample,y=Reads,fill=Bacteria))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
