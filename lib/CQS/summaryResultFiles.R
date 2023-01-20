options(bitmapType='cairo')

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
canBeEmpty<-strsplit(fileList$CanFileEmpty, ",")
fileSizeRaw<-lapply(temp,function(x) file.size1(x))
fileSize<-lapply(fileSizeRaw,addUnitToSize)
fileSizeTotalRaw<-sapply(fileSizeRaw,sum,na.rm=T)
fileSizeTotal<-addUnitToSize(fileSizeTotalRaw)

fileExistPercent<-sapply(c(1:length(fileSizeRaw)),function(i){
  df<-data.frame(Size=fileSizeRaw[[i]], Empty=canBeEmpty[[i]])
  df<-df[!(is.na(df$Size) & df$Empty=="True"),]
  df$RequiredSize<-lapply(df$Empty, function(x){
    if(x == "True"){
      return(0)
    }else{
      return(100)
    }
  })
  
  perc<-length(which(df$Size>=df$RequiredSize&!is.na(df$Size)))/length(df$Size)
  return(perc)
})

Result<-rep("FAIL",length(fileExistPercent))
Result[which(fileExistPercent==1)]<-"PASS"
Result[which(fileExistPercent>0 & fileExistPercent<1)]<-"WARN"
fileSize<-sapply(fileSize,function(x) paste(x,collapse=","))

ResultOut<-data.frame(fileList,FileSize=fileSize,FileSizeTotalRaw=fileSizeTotalRaw,FileSizeTotal=fileSizeTotal,FileExistPercent=fileExistPercent,Result=Result,stringsAsFactors=FALSE)
write.csv(ResultOut,paste0(fileListName,".check.csv"))

colors<-c("light green", "skyblue", "red")
names(colors)<-c("PASS","WARN","FAIL")

for (step in unique(ResultOut$StepName)) {
  tableForPlot<-ResultOut[which(ResultOut$StepName==step),]
  #remove some smallRNA tasks from scatter
  tableForPlot<-tableForPlot[!grepl('_bacteria.\\d{3}', tableForPlot$SampleName),]

  failLength<-length(grep("FAIL",tableForPlot$Result))
  warnLength<-length(grep("WARN",tableForPlot$Result))
  if (failLength>0) {
    print(paste0("There are ",failLength," FAIL in ",step,"."))
    message<-as.vector(apply(tableForPlot[grep("FAIL",tableForPlot$Result),c("TaskName","SampleName")],1,function(x) paste(x,collapse=": ")))
    if (length(message)<=10) {
      print(message)
    } else {
      print("First 10 tasks:")
      print(message[1:10])
    }
  }
  if (warnLength>0) {
    print(paste0("There are ",warnLength," WARN in ",step,"."))
    message<-as.vector(apply(tableForPlot[grep("WARN",tableForPlot$Result),c("TaskName","SampleName")],1,function(x) paste(x,collapse=": ")))
    if (length(message)<=10) {
      print(message)
    } else {
      print("First 10 tasks:")
      print(message[1:10])
    }
  }
  if (failLength==0 & warnLength==0) {
    print(paste0("All tasks are successfully finished in ",step,"."))
  }
  tableForPlot$TaskName<-factor(tableForPlot$TaskName,levels=rev(unique(tableForPlot$TaskName)))
  tableForPlot$Result<-factor(tableForPlot$Result,levels=c("PASS","WARN","FAIL"))
  width=max(3000, 80 * length(unique(tableForPlot$SampleName)))
  height=max(2000, 100 * length(unique(tableForPlot$TaskName)))
  png(file=paste0(fileListName,"_",step,".png"), height=height, width=width, res=300)
  g<-ggplot(tableForPlot, aes(SampleName, TaskName))+
      geom_tile(data=tableForPlot, aes(fill=Result), color="white") +
      scale_fill_manual(values=colors) +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold"),
          axis.text.y = element_text(size=11, face="bold")) +
      labs(title = paste0(step,": Status"))+
      coord_equal()
  print(g)
  dev.off()
  
  temp<-aggregate(tableForPlot$FileSizeTotalRaw, list(factor(tableForPlot$TaskName)), function(x) if (median(x,na.rm=T)==0) {return(max(x,na.rm=T))} else {median(x,na.rm=T)})
  taskFileSizeMedian<-temp[,2]
  names(taskFileSizeMedian)<-temp[,1]
  
  tableForPlot$Task<-paste0(tableForPlot$TaskName," (",addUnitToSize(taskFileSizeMedian)[tableForPlot$TaskName],")")
  tableForPlot$Task<-factor(tableForPlot$Task,levels=rev(unique(tableForPlot$Task)))
  tableForPlot$Log2RelativeSize<-log2(tableForPlot$FileSizeTotalRaw/taskFileSizeMedian[tableForPlot$TaskName])
  
  if (all(is.na(tableForPlot$Log2RelativeSize))) {
    print(paste0("All tasks in ",step," have 0 file size."))
  } else {
    png(file=paste0(fileListName,"_",step,".RelativeFileSize.png"),height=height, width=width, res=300)
    g<-ggplot(tableForPlot, aes(SampleName, Task))+
        geom_tile(data=tableForPlot, aes(fill=Log2RelativeSize), color="white") +
        scale_fill_gradient2(low="lightgreen", high="red") +
        theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0.5, face="bold"),
            axis.text.y = element_text(size=11, face="bold")) +
        labs(title = paste0(step,": Relative File Size"))+
        coord_equal()
    print(g)
    dev.off()
  }

}


