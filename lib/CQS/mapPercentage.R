options(bitmapType='cairo')

# TODO: Add comment
# 
# Author: Quanhu Sheng
###############################################################################

resultFile<-outFile
readFileList<-parSampleFile1

require(XML)
library(ggplot2)

readFiles<-read.delim(readFileList,header=F,as.is=T)
samples<-unique(readFiles$V2)
curSample<-samples[1]
for (curSample in samples) {
  resFile<-paste0(curSample, ".csv")
  if(file.exists(resFile)){
    df<-read.csv(resFile,check.names=F)
  }else{
    curFiles<-readFiles[readFiles$V2==curSample, ]
    curFiles$Category<-"Nonhost"
    curFiles$Category[2]<-"Host smallRNA"
    curFiles$Category[3]<-"Host genome"
    allReads<-read.table(curFiles[1,1], sep="\t", header=T, row.names=1)[,c("Count"),drop=F]
    allReads$Mapped<-"Unmapped"
    mapIndex<-3
    for (mapIndex in c(2:nrow(curFiles))){
      mapFile<-curFiles[mapIndex,1]
      if(grepl(".xml", mapFile)){
        cat("Xml: ", mapFile, "\n")
        xml_data <- xmlToList(mapFile)
        temp<-xml_data[["queries"]]
        qnames<-sapply(temp,function(x) {
          xLength=length(x);
          x[[xLength]]["name"];
        })
        qnames<-gsub(":CLIP_.*", "", qnames )
      }else{
        cat("Txt: ", mapFile, "\n")
        qnames<-read.table(mapFile, sep="\t", header=T, stringsAsFactor=F)$Query
      }
      allReads[qnames,"Mapped"]<-curFiles[mapIndex,"Category"]
    }
    
    uniqueCounts<-unique(allReads$Count)
    x<-1
    res<-lapply(uniqueCounts, function(x){
      ta<-table(allReads$Mapped[allReads$Count==x])
      df<-data.frame(ta)
      df$Count=x
      return(df)
    })
    df<-do.call(rbind, res)
    colnames(df)<-c("Category", "Value", "ReadCount")
    
    write.csv(df, paste0(curSample, ".csv"))
  }
  
  df<-df[df$ReadCount <=20,]
  df$ReadCount<-factor(df$ReadCount)
  df$Measure<-"Frequency"
  
  df2<-df
  temp<-tapply(df[,"Value"],df[,"ReadCount"],sum)
  df2$Value<-df[,"Value"]/temp[df[,"ReadCount"]] * 100
  df2$Measure<-"Percentage"
  
  dfall<-rbind(df, df2)
  dfall$Category=factor(dfall$Category, levels=c("Host smallRNA","Host genome","Nonhost","Unmapped"))
  g<-ggplot(dfall, aes(x=ReadCount, y=Value, fill=Category)) + geom_bar(stat="identity") + facet_wrap(~Measure, nrow=2, scales="free_y")
  png(file=paste0(curSample, ".mapped.png"), width=1600, height=2000, res=300)
  print(g)
  dev.off()
}
