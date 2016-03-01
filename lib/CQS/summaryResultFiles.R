fileListName<-parFile1

#fileListName<-commandArgs()[7]
#fileListName<-"/scratch/cqs/zhaos/vickers/20160123_smallRNA_3018-KCV-55_human/sequencetask/result/3018-KCV-55_st_expect_result.tsv"
#setwd("/scratch/cqs/zhaos/vickers/20160123_smallRNA_3018-KCV-55_human/sequencetask/result")

library(ggplot2)
file.size1<-function(x) {
	starIndex<-grep("\\*",x)
	if (length(starIndex)>0) {
		result<-rep(NA,length(x))
		for (i in starIndex) {
			basenameStar<-basename(x[i])
			dirnameStar<-dirname(x[i])
			temp<-list.files(dirnameStar,basenameStar,full.names=T)
			result[i]<-sum(file.size(temp))
		}
		result[-starIndex]<-file.size(x[-starIndex])
		return(result)
	} else {
		return(file.size(x))
	}
}
addUnitToSize<-function(x) {
	xCut<-cut(x,c(-Inf,1024,1024^2,1024^3,Inf))
	reSizeFactor<-c(1,1024,1024^2,1024^3)
	reSizeUnit<-c("B","K","M","G")
	result<-x/reSizeFactor[as.integer(xCut)]
	result<-paste0(round(result,0),reSizeUnit[as.integer(xCut)])
	result[is.na(x)]<-NA
	
	return(result)
}

fileList<-read.delim(fileListName,as.is=T,header=T)
temp<-strsplit(fileList$FileList,",")
fileSizeRaw<-lapply(temp,function(x) file.size1(x))
fileSize<-lapply(fileSizeRaw,addUnitToSize)
fileSizeTotalRaw<-sapply(fileSizeRaw,sum,na.rm=T)
fileSizeTotal<-addUnitToSize(fileSizeTotalRaw)
fileExistPercent<-sapply(fileSizeRaw,function(x) length(which(x>100&!is.na(x)))/length(x))
Result<-rep("FAIL",length(fileExistPercent))
Result[which(fileExistPercent==1)]<-"PASS"
Result[which(fileExistPercent>0 & fileExistPercent<1)]<-"WARN"
fileSize<-sapply(fileSize,function(x) paste(x,collapse=","))

ResultOut<-data.frame(fileList,FileSize=fileSize,FileSizeTotalRaw=fileSizeTotalRaw,FileSizeTotal=fileSizeTotal,FileExistPercent=fileExistPercent,Result=Result,stringsAsFactors=FALSE)
write.csv(ResultOut,paste0(fileListName,".check.csv"))

for (step in unique(ResultOut$StepName)) {
	tableForPlot<-ResultOut[which(ResultOut$StepName==step),]
	failLength<-length(grep("FAIL",tableForPlot$Result))
	warnLength<-length(grep("WARN",tableForPlot$Result))
	if (failLength>0) {
		print(paste0("There are ",failLength," FAIL in ",step,"."))
	}
	if (warnLength>0) {
		print(paste0("There are ",warnLength," WARN in ",step,"."))
	}
	if (failLength==0 & warnLength==0) {
		print(paste0("All tasks are successfully finished in ",step,"."))
	}
	tableForPlot$TaskName<-factor(tableForPlot$TaskName,levels=rev(unique(tableForPlot$TaskName)))
	tableForPlot$Result<-factor(tableForPlot$Result,levels=c("PASS","WARN","FAIL"))
	width=max(2500, 60 * length(unique(tableForPlot$SampleName)))
	height=max(2000, 60 * length(unique(tableForPlot$TaskName)))
	png(file=paste0(fileListName,"_",step,".png"), height=height, width=width, res=300)
	g<-ggplot(tableForPlot, aes(SampleName, TaskName))+
			geom_tile(data=tableForPlot, aes(fill=Result), color="white") +
			scale_fill_manual(values=c("light green", "skyblue", "red")) +
			theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold"),
					axis.text.y = element_text(size=11, face="bold")) +
			coord_equal()
	print(g)
	dev.off()
	
	temp<-aggregate(tableForPlot$FileSizeTotalRaw, list(factor(tableForPlot$TaskName)), median,na.rm=T)
	taskFileSizeMedian<-temp[,2]
	names(taskFileSizeMedian)<-temp[,1]
	
	tableForPlot$Task<-paste0(tableForPlot$TaskName," (",addUnitToSize(taskFileSizeMedian)[tableForPlot$TaskName],")")
	tableForPlot$RelativeSize<-tableForPlot$FileSizeTotalRaw/taskFileSizeMedian[tableForPlot$TaskName]
	
	png(file=paste0(fileListName,"_",step,".RelativeFileSize.png"),height=height, width=width, res=300)
	g<-ggplot(tableForPlot, aes(SampleName, Task))+
			geom_tile(data=tableForPlot, aes(fill=RelativeSize), color="white") +
			theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0.5, face="bold"),
					axis.text.y = element_text(size=11, face="bold")) +
			coord_equal()
	print(g)
	dev.off()
}


