resultFile<-outFile
countFiles<-parSampleFile1

library(reshape)
library(ggplot2)

resultTable<-NULL
for (countFile in countFiles[1,]) {
	if (grepl(".csv$",countFile)) {
		countTable<-read.csv(countFile,header=T,row.names=1)
	} else {
		countTable<-read.delim(countFile,header=T,row.names=1)
	}
	resultTable<-rbind(resultTable,colSums(countTable))
}


row.names(resultTable)<-countFiles[2,]
dataForFigure<-melt(resultTable)
colnames(dataForFigure)<-c("Task","Sample","Reads")
write.csv(dataForFigure,paste0(resultFile,".TaskReads.csv"))

ggplot(dataForFigure,aes(x=Sample,y=Reads))+geom_bar(stat="identity",position="dodge",aes(fill=Task))

