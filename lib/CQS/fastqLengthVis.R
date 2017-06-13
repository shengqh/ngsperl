options(bitmapType='cairo')

resultFile<-outFile
lengthFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3

if(!exists("visLayoutAlphabet")){
  visLayoutAlphabet<-FALSE
}

library(reshape2)
library(ggplot2)
library(cowplot)

#source("/home/zhaos/source/r_cqs/vickers/codesToPipeline/countTableVisFunctions.R")

lengthFiles<-read.delim(lengthFileList,as.is=T,header=F)

lengthAll<-NULL
for (i in 1:nrow(lengthFiles)) {
  lengthFile<-lengthFiles[i,1]
  sampleName<-lengthFiles[i,2]
  if(exists("ignoreSamples") && sampleName %in% ignoreSamples){
    cat("Sample ", sampleName, " ignored.\n")
    next
  }
  lengthOne<-read.delim(lengthFile,as.is=T,header=T)
  lengthOne$Sample<-sampleName
  lengthAll<-rbind(lengthAll,lengthOne)
}

lengthAllTable<-acast(lengthAll,Len~Sample,value.var="Count")
write.csv(lengthAllTable,paste0(resultFile,".csv"), row.names=F)

cellWidth=800
scales="free_x"
if(exists("free_y") && free_y){
  cellWidth=1000
  scales="free"
}

width=max(2000, cellWidth * (ceiling(sqrt(ncol(lengthAllTable)))+1))
png(paste0(resultFile,".png"),width=width,height=width,res=300)
p=ggplot(lengthAll, aes(x=Len, y=Count)) + 
  geom_bar(stat="identity", width=.5) + 
  facet_wrap(~Sample, scales=scales) +
  xlab("Read length") + 
  ylab("Read count")
print(p)
dev.off()

if (groupFileList!="") {
  sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
  #keep the groups with samples in the count table
  sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(lengthAllTable)),]
  
  datBySampleGroup<-mergeTableBySampleGroup(lengthAllTable,sampleToGroup)
  temp<-melt(datBySampleGroup)
  colnames(temp)<-c("Position","Group","Percent")
  
  if (groupVisLayoutFileList != ""){
    visLayout<-readLayout(groupVisLayoutFileList,visLayoutAlphabet)
    temp<-data.frame(temp,visLayout[temp$Group,])
    width<-getLayoutWidth(visLayout, 2000, cellWidth)
    png(paste0(resultFile,".group.png"),width=width,height=width,res=300)
    p<-ggplot(temp,aes(x=Position,y=Percent))+geom_line()
    p<-p+facet_grid(Row_Group~Col_Group, scales="free")
  }
  else{
    p<-ggplot(temp,aes(x=Position,y=Percent,colour= Group))+geom_line()
    
    if(exists("facet") && facet){
      width=max(2000, cellWidth * (ceiling(sqrt(ncol(temp)))+1))
      png(paste0(resultFile,".group.png"),width=width,height=width,res=300)
      p<-p+facet_grid(Row_Group~Col_Group, scales="free")
    }else{
      png(paste0(resultFile,".group.png"),width=2000,height=1500,res=300)
    }  
  }
  print(p)
  dev.off()
}
