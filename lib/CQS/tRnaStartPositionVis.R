# TODO: Add comment
# 
# Author: shili
###############################################################################


options(bitmapType='cairo')
maxFeature=40

groupFileList=parSampleFile1
visLayoutFileList=parSampleFile2
positionFile = parFile1
totalCountFile<-parFile2


#load Rcpp package first because of the error with reshape2 package
library(Rcpp)
library(grid)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

if (!exists("countName")){
	countName = "MappedReads"
}


tRnaName2Group<-function(x,toAnticodn=TRUE) {
	x=gsub("^tRNA:","",x)
	x=gsub("^nm-|^nmt-","",x)
	
	temp=strsplit(x,"-")
	if (toAnticodn) {
		tRnaGroup=sapply(temp,function(x) paste(x[2:3],collapse=""))
	} else {
		tRnaGroup=sapply(temp,function(x) x[2])
	}
	return(tRnaGroup)
}


position<-read.delim(positionFile, header=T,as.is=T)
samples<-unique(position$File)

groupInfo<-read.delim(groupFileList, header=F, as.is=T)
groups<-unique(groupInfo$V2)
groupSize<-table(groupInfo[,2])

totalCount<-read.delim(totalCountFile,as.is=T,header=T,row.names=1,check.names=FALSE)
totalCount<-unlist(totalCount[countName,])

position$FeatureLabel<-paste0(position$Feature,"(",round(position$TotalCount,0),")")
position$PositionCountFraction<-as.vector(position$PositionCount/totalCount[position$File])

allPosition<-NULL
for(group in groups){
	groupSamples<-groupInfo$V1[groupInfo$V2==group]
	groupPositions<-position[position$File %in% groupSamples,]
	groupPositions$Group<-group
	allPosition<-rbind(allPosition, groupPositions)
}

#tRNA Group
allPosition$smallRNAGroup<-tRnaName2Group(allPosition$Feature,toAnticodn=TRUE)
smallRNAGroupSize<-tapply(allPosition$Feature,allPosition$smallRNAGroup,function(x) length(unique(x)))

if (length(unique(allPosition$Feature))>maxFeature) {
	allPositionUniqueFeature<-unique(allPosition[,c("File","Feature","TotalCount")])
	featureTotalCount<-tapply(allPositionUniqueFeature$TotalCount,allPositionUniqueFeature$Feature,sum)
	featureKept<-names(featureTotalCount)[rev(order(featureTotalCount))][1:maxFeature]
	allPositionKept<-allPosition[which(allPosition$Feature %in% featureKept),]
	allPosition<-allPositionKept
}


allPositionByGroup<-aggregate(x = allPosition, by = list(allPosition$Feature,allPosition$Group, allPosition$Position), FUN = function(x) if(is.numeric(x)| is.integer(x)) {sum(x)} else {x[1]})
allPositionByGroup$GroupPositionCount<-as.vector(allPositionByGroup$PositionCount/groupSize[allPositionByGroup$Group])
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
} else {
	allPositionByGroup$Group<-factor(allPositionByGroup$Group,levels=groups)
}

featureNumber<-length(unique(allPositionByGroup$Feature))

xRange<-range(allPositionByGroup$Position)
if (xRange[1]>=-5) {
	xRange[1]=-5
}
if (xRange[2]>=110) {
	xRange[2]=110
}


axisTextSize=16
stripTextSize=12
height=max(featureNumber*200,3000)
groupCount=length(unique(allPositionByGroup$Group))
width=max(groupCount*1000, 3000)

png(paste0(outFile,".barplot.png"),width=width,height=height,res=300)
m <- ggplot(allPositionByGroup, aes(x = Position,y=GroupPercentage)) +
		geom_bar(stat="identity") +
		theme_bw()+
		ylab("Read fraction (read counts/total reads)")+
		theme(text = element_text(size=axisTextSize))+
		theme(legend.text = element_text(size=axisTextSize))+
		theme(strip.text.y = element_text(size=stripTextSize,angle = 0))
#		scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(featureNumber)) + 
#		xlim(-10, maxPos+5) +
#		theme(legend.key.size = unit(0.4, "cm"), legend.position="right") +
#		guides(fill= guide_legend(ncol=ncols,keywidth=1, keyheight=1.5))

#m <- ggplot(allPositionByGroup, aes(x = Position,y=GroupPositionCountFraction))+geom_bar(stat="identity")
print(m+facet_grid(smallRNAGroup~Group,scale="free"))
dev.off()



