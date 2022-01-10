source("countTableVisFunctions.R")
options(bitmapType='cairo')

resultFile<-outFile
lengthFileList<-parSampleFile1
groupFileList<-parSampleFile2
groupVisLayoutFileList<-parSampleFile3
free_y=1

if(!exists("visLayoutAlphabet")){
  visLayoutAlphabet<-FALSE
}

library(reshape2)
library(ggplot2)

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
write.csv(lengthAllTable,paste0(resultFile,".csv"), row.names=T)

cellWidth=1500
scales="free"

textTitle<-element_text(face= "bold", color = "black", size=22, hjust=0.5)
text20Bold<-element_text(face= "bold", color = "black", size=20)
text20<-element_text(color = "black", size=20)

sampleCount<-ncol(lengthAllTable)
facetColCount = ceiling(sqrt(sampleCount))
facetRowCount = ceiling( sampleCount * 1.0 / facetColCount)
width=max(2000, cellWidth * facetColCount) + 100
height=max(2000, cellWidth * facetRowCount)

png(paste0(resultFile,".png"),width=width,height=height,res=300)

p=ggplot(lengthAll, aes(x=Len, y=Count)) + 
  geom_bar(stat="identity", width=.5) + 
  theme_bw3() + 
  labs(x = "Read Length", 
       y = "Read Count")+
  theme(plot.title = textTitle,
        axis.title = text20Bold,
        axis.text = text20,
        axis.line = element_line(colour = "gray75", size=0.73, linetype = "solid"),
        axis.ticks = element_line(size=0.73),axis.ticks.length=unit(0.3,"cm"),
        strip.text = text20Bold,
        legend.text = text20Bold,
        legend.title = element_blank()) +
        facet_wrap(~Sample, ncol=facetColCount, scales=scales)
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
  p<-p+theme_bw3()
  print(p)
  dev.off()
}
