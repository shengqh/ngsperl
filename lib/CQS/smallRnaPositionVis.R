# TODO: Add comment
# 
# Author: zhaos
###############################################################################

maxFeature=25

groupFileList=parSampleFile1
visLayoutFileList=parSampleFile2
positionFile = parFile1
totalCountFile<-parFile2

if(!exists("useLeastGroups")){
  useLeastGroups<-FALSE
}

#load Rcpp package first because of the error with reshape2 package
library(Rcpp)

library(reshape2)
library(ggplot2)
library(RColorBrewer)


smallRNAGrouping<-function(x) {
  tRNAHeaders<-c("Glu","Gly","Asp","Val","His","Thr","Leu","Lys","Cys","Pro",
      "Tyr","Ala","SeC","Asn","Ser","Trp","Gln","iMet","Arg","Ile",
      "Met","Phe","Sup")
  tRNAHeadersPattern<-paste0(tRNAHeaders,collapse="|")
  if (all(grepl("^RNU|^RNVU|^U\\d+",head(as.character(x))))) {
    return(1) #group-able snRNA
  } else if (all(head(as.character(x)) %in% tRNAHeaders)) {
    return(2) #tRNA headers only (Gly)
  } else if (all(grepl(tRNAHeadersPattern,head(as.character(x))))) {
    return(3) #tRNA with headers and anticodn only (GlyACC)
  }else { #NOT group-able snRNA
    return(0)
  }
}

smallRnaName2Group<-function(x,groupSnRNA=1) {
  if (groupSnRNA==1) { #snRNA
    snRnaGroup<-sapply(strsplit(as.character(x),"-"),function(y) y[1])
    snRnaGroup<-gsub("RNVU","RNU",snRnaGroup)
    snRnaGroup<-gsub("([A-Z]+[0-9]+)[A-Z]+","\\1",snRnaGroup)
  } else if (groupSnRNA==3) { #tRNA with headers and anticodn only (GlyACC)
    snRnaGroup<-substr(x,0,3)
  } else {
    snRnaGroup<-x
  }
  return(snRnaGroup)
}

allPosition<-read.delim(positionFile, header=T,as.is=T)
samples<-unique(allPosition$File)

groupInfo<-getSampleInGroup(groupFileList, samples, useLeastGroups=useLeastGroups)
groupSize<-table(groupInfo[,1])
totalCount<-read.delim(totalCountFile,as.is=T,header=T,row.names=1)
totalCount<-unlist(totalCount["MappedReads",])

allPosition$FeatureLabel<-paste0(allPosition$Feature,"(",round(allPosition$TotalCount,0),")")
allPosition$Group<-groupInfo[allPosition$File,1]
allPosition$PositionCountFraction<-as.vector(allPosition$PositionCount/totalCount[allPosition$File])

doSmallRNAGrouping<-smallRNAGrouping(unique(allPosition$Feature))
if (doSmallRNAGrouping==1 | doSmallRNAGrouping==3) {
  allPosition$smallRNAGroup<-smallRnaName2Group(allPosition$Feature,doSmallRNAGrouping)
}

if (length(unique(allPosition$Feature))>maxFeature) {
  allPositionUniqueFeature<-unique(allPosition[,c("File","Feature","TotalCount")])
  featureTotalCount<-tapply(allPositionUniqueFeature$TotalCount,allPositionUniqueFeature$Feature,sum)
  featureKept<-names(featureTotalCount)[rev(order(featureTotalCount))][1:maxFeature]
  allPositionKept<-allPosition[which(allPosition$Feature %in% featureKept),]
  allPosition<-allPositionKept
}

allPositionByGroup<-aggregate(x = allPosition, by = list(allPosition$Feature,allPosition$Group, allPosition$Position), FUN = function(x) if(is.numeric(x)| is.integer(x)) {sum(x)} else {x[1]})
allPositionByGroup$GroupPercentage<-as.vector(allPositionByGroup$Percentage/groupSize[allPositionByGroup$Group])
allPositionByGroup$GroupPositionCountFraction<-as.vector(allPositionByGroup$PositionCountFraction/groupSize[allPositionByGroup$Group])
allPositionByGroup$Position<-allPositionByGroup$Group.3

if (visLayoutFileList!="") {
  visLayout<-read.delim(visLayoutFileList,as.is=T,header=F)
  visLayout<-sapply(split(visLayout[,1],visLayout[,2]),function(x) x)
  row.names(visLayout)<-visLayout[,"Groups"]
  visLayout<-data.frame(visLayout[,-which(colnames(visLayout)=="Groups")])
  visLayout$Col_Group<-factor(visLayout$Col_Group,levels=unique(visLayout$Col_Group))
  visLayout$Row_Group<-factor(visLayout$Row_Group,levels=unique(visLayout$Row_Group))
  
  allPositionByGroup<-data.frame(allPositionByGroup,visLayout[allPositionByGroup[,"Group"],])
}
featureNumber<-length(unique(allPositionByGroup$Feature))

if (visLayoutFileList!="") {
  height=max(length(unique(allPositionByGroup$Row_Group))*featureNumber*80,3000)
  width=max(length(unique(allPositionByGroup$Col_Group))*2000,3000)
} else {
  height=max(featureNumber*100,3000)
  width=max(length(unique(allPositionByGroup$Group))*2000,3000)
}
xRange<-range(allPositionByGroup$Position)
if (xRange[1]>=-5) {
  xRange[1]=-10
} else {
  xRange[1]=xRange[1]-5
}
if (xRange[2]>=110) {
  xRange[2]=120
} else {
  xRange[2]=xRange[2]+10
}

axisTextSize=12
stripTextSize=12
png(paste0(outFile,".png"), width=width, height=height, res=300)
p<-ggplot(allPositionByGroup,aes(x=Position,y=Feature,size=GroupPercentage,colour=GroupPercentage))+
    geom_point()+
    scale_size_continuous(range = c(0.1,3))+
    scale_colour_gradient(low="indianred1",high="darkred")+
    xlim(xRange)+ 
    theme_bw()+
    theme(text = element_text(size=axisTextSize),axis.text = element_text(size=axisTextSize),
        axis.title = element_text(size=axisTextSize),
        strip.text.x = element_text(size=stripTextSize),
        strip.text.y = element_text(size=stripTextSize,angle = 0))+
    theme(legend.position="none")
if (visLayoutFileList!="" & "smallRNAGroup" %in% colnames(allPositionByGroup)) {
  p<-p+facet_grid(Row_Group+smallRNAGroup~Col_Group,space = "free",scale="free")
} else if (visLayoutFileList!="") {
  p<-p+facet_grid(Row_Group~Col_Group,space = "free",scale="free")
} else {
  p<-p+facet_grid(.~Group,space = "free",scale="free")
}
print(p)
dev.off()

if (visLayoutFileList!="") {
  height=max(length(unique(allPositionByGroup$Row_Group))*2000,3000)
  width=max(length(unique(allPositionByGroup$Col_Group))*2000,3000)
} else {
  height=3000
  width=max(length(unique(allPositionByGroup$Group))*2000,3000)
}
maxPos<-max(allPositionByGroup$Position)
png(paste0(outFile,".allPositionBar.png"),width=width,height=height,res=300)
m <- ggplot(allPositionByGroup, aes(x = Position,y=GroupPositionCountFraction,fill=Feature)) +
    geom_bar(stat="identity") +
    theme_bw()+
    theme(legend.key.size = unit(0.4, "cm"))+
    ylab("cumulative read fraction (read counts/total reads)")+
    theme(text = element_text(size=20))+theme(legend.text = element_text(size=16))+
    guides(fill= guide_legend(ncol = 1,keywidth=1, keyheight=1.5))+
    scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(featureNumber)) + 
    xlim(-10, maxPos+5)
if (visLayoutFileList!="") {
  m<-m+facet_grid(Row_Group~Col_Group,space = "free",scale="free")
} else {
  m<-m+facet_grid(.~Group,space = "free",scale="free")
}
print(m)
dev.off()



