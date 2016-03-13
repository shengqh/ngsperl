resultFile<-outFile
categoryFileList<-parSampleFile1

library(reshape2)


categoryFiles<-read.delim(categoryFileList,as.is=T,header=F)

categoryAll<-NULL
categoryFigure<-NULL
for (i in 1:nrow(categoryFiles)) {
	categoryFile<-categoryFiles[i,1]
	sampleName<-categoryFiles[i,2]
	categoryOne<-read.delim(categoryFile,as.is=T,header=F,row.names=1,comment.char = "#")
	categoryOne$Sample<-sampleName
	
	Unmapped<-categoryOne["TotalReads",1]-categoryOne["MappedReads",1]
	smallRna<-sum(categoryOne[-c(1:3),1])
	otherMapped<-categoryOne["MappedReads",1]-smallRna
	categoryOneFigure<-data.frame(Category=c("Unmapped","Other Mapped","Small RNA"),Count=c(Unmapped,otherMapped,smallRna),Sample=sampleName)
	
	categoryAll<-rbind(categoryAll,categoryOne)
	categoryFigure<-rbind(categoryFigure,categoryOneFigure)
}

png(paste0(resultFile,".png"), width=4000, height=3000, res=300)
g<-ggplot(gcounts,aes(x=Sample, y=Count,fill=Category))+geom_bar(stat="identity")+
		theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
		theme(legend.position = "top")
print(g)

dev.off()











