###############################################################################
# Funtions in count table barplot and pie chart
###############################################################################
library(reshape2)
library(ggplot2)
library(RColorBrewer)

orderDataByNames<-function(x,orderKey,orderNames) {
	orderedInd<-NULL
	for (i in orderNames) {
		orderedInd<-c(orderedInd,which(orderKey==i))
	}
	return(x[orderedInd,])
}

mergeTableBySampleGroup<-function(x,sampleToGroup) {
	xRatio<-t(t(x)/colSums(x))
	groupLength<-length(unique(sampleToGroup[,2]))
	xRatioGroupMean<-matrix(NA,ncol=groupLength,nrow=nrow(x))
	colnames(xRatioGroupMean)<-unique(sampleToGroup[,2])
	row.names(xRatioGroupMean)<-row.names(x)
	for (i in 1:groupLength) {
		currentSample<-sampleToGroup[which(sampleToGroup[,2]==colnames(xRatioGroupMean)[i]),1]
		xRatioGroupMean[,i]<-rowMeans(xRatio[,currentSample])
	}
	return(xRatioGroupMean)
}

basicPie<-function(x,maxCategory=10,main="",addPercent=F) {
	nameToCountForFigure<-na.omit(rev(sort(x)))
	if (length(nameToCountForFigure)>maxCategory) {
		nameToCountForFigure<-c(nameToCountForFigure[1:(maxCategory-1)],Other=sum(nameToCountForFigure[-(1:(maxCategory-1))]))
	}
	if (addPercent) {
		temp<-round(nameToCountForFigure/sum(nameToCountForFigure)*100,0)
		names(nameToCountForFigure)<-paste0(names(nameToCountForFigure)," ",temp,"%")
	}
	pie(nameToCountForFigure,col=rainbow(length(nameToCountForFigure)),main=main)
}

makeColors<-function(n,colorNames="Set1") {
	maxN<-brewer.pal.info[colorNames,"maxcolors"]
	if (n<=maxN) {
		colors<-brewer.pal(n, colorNames)
	} else {
		colors<-colorRampPalette(brewer.pal(maxN, colorNames))(n)
	}
	return(colors)
}

tableMaxCategory<-function(dat,maxCategory=NA) {
	if (!is.na(maxCategory) & nrow(dat)>maxCategory) {
		temp<-apply(dat,2,function(x) rev(order(x))[1:maxCategory])
		categoryKeptInd<-sort(unique(as.vector(temp)))
		datForFigure<-dat[categoryKeptInd,]
		if (length(categoryKeptInd)<nrow(dat)) {
			datForFigure<-rbind(datForFigure,Other=colSums(dat[-categoryKeptInd,,drop=FALSE]))
		}
	} else {
		datForFigure<-dat
	}
	return(datForFigure)
}

tableBarplot<-function(dat,maxCategory=5,x="Sample", y="Reads",fill="Category",facet=NA,varName=if (is.na(facet)) c(fill,x,y) else c(facet,x,y),transformTable=TRUE,textSize=20,ylab=y,colorNames="Set1") {
	if (transformTable) {
		datForFigure<-tableMaxCategory(dat,maxCategory=maxCategory)
		
#		datForFigure$Groups<-row.names(dat)
		datForFigure<-melt(as.matrix(datForFigure))
		colnames(datForFigure)<-varName
	} else {
		datForFigure<-dat
	}
	if (!is.na(fill)) {
		p<-ggplot(datForFigure,aes_string(x=x,y=y,fill=fill))
		if (length(unique(datForFigure[,fill]))<=7 & sum(nchar(as.character(unique(datForFigure[,fill]))))<70) {
			p<-p+theme(legend.position = "top")+
					guides(fill = guide_legend(nrow = 1,keywidth = 2, keyheight = 2))
		} else {
			p<-p+	guides(fill = guide_legend(ncol = 1,keywidth = 1, keyheight = 1))
		}
		if (colorNames!="") {
			colors<-makeColors(length(unique(datForFigure[,fill])),colorNames)
			p<-p+scale_fill_manual(values=colors)
		}
	} else if (!is.na(facet)) {
		p<-ggplot(datForFigure,aes_string(x=x,y=y))+facet_wrap(c(facet))
	} else {
		p<-ggplot(datForFigure,aes_string(x=x,y=y))
	}
	p<-p+geom_bar(stat="identity")+
#			guides(fill= guide_legend(title = groupName))+
			theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
			theme(axis.text = element_text(size=textSize),legend.text=element_text(size=textSize),
					axis.title = element_text(size=textSize),legend.title= element_text(size=textSize))+
			ylab(ylab)

	return(p)
}

tableBarplotToFile<-function(dat,fileName,totalCountFile="",maxCategory=5,textSize=9,transformTable=T,height=1500,...) {
	if (totalCountFile!="") { #normlize with total count *10^6
		totalCount<-read.csv(totalCountFile,header=T,as.is=T,row.names=1,check.names=FALSE)
		totalCount<-unlist(totalCount["Reads for Mapping",])
		dat<-10^6*t(t(dat)/totalCount[colnames(dat)])
		ylab<-"Mapped Reads per Million"
	} else {
		ylab<-"Reads"
	}
	width<-max(3000,75*ncol(dat))
	height<-height
	png(fileName,width=width,height=height,res=300)
	p<-tableBarplot(dat,maxCategory=maxCategory,textSize=textSize,ylab=ylab,transformTable=transformTable,...)
	print(p)
	dev.off()
}

#changed from function in http://mathematicalcoffee.blogspot.com/2014/06/ggpie-pie-graphs-in-ggplot2.html
ggpie <- function (dat, fill="Category", y="Reads",facet="Sample", maxCategory=NA,main=NA, percent=T,textSize=15,colorNames="Set1",transformTable=TRUE,reOrder=TRUE) {
	if (transformTable) {
		datForFigure<-tableMaxCategory(dat,maxCategory=maxCategory)
		
	    if (percent) {
		    datForFigure<-t(t(datForFigure)/colSums(datForFigure))
	    }
		if (reOrder) {
			if (row.names(datForFigure)[nrow(datForFigure)]=="Other") {
				categoryOrderedNames<-c(row.names(datForFigure)[-nrow(datForFigure)][rev(order(rowSums(datForFigure[-nrow(datForFigure),])))],"Other")
			} else {
				categoryOrderedNames<-row.names(datForFigure)[rev(order(rowSums(datForFigure)))]
			}
		} else {
			categoryOrderedNames<-row.names(datForFigure)
		}
		datForFigure<-melt(as.matrix(datForFigure),as.is=T)
		colnames(datForFigure)<-c(fill,facet,y)
		datForFigure[,fill]<-factor(datForFigure[,fill],levels=categoryOrderedNames)
		datForFigure<-orderDataByNames(datForFigure,datForFigure[,fill],categoryOrderedNames)
	} else {
		datForFigure<-dat
		if (percent) {
			if (!is.na(facet)) {
				temp<-tapply(datForFigure[,y],datForFigure[,facet],sum)
				datForFigure[,y]<-datForFigure[,y]/temp[datForFigure[,facet]]
			} else {
				datForFigure[,y]<-datForFigure[,y]/sum(datForFigure[,y])
			}
		}
		categoryOrderedNames<-levels(datForFigure[,fill])[rev(order(tapply(datForFigure[,y],datForFigure[,fill],sum)))]
		datForFigure[,fill]<-factor(datForFigure[,fill],levels=categoryOrderedNames)
		datForFigure<-orderDataByNames(datForFigure,datForFigure[,fill],categoryOrderedNames)
	}
	p = ggplot(datForFigure, aes_string(x=factor(1), y=y, fill=fill)) +
			geom_bar(width=1, stat='identity', color='black') +
			guides(fill=guide_legend(keywidth = 1.5, keyheight = 1.5,override.aes=list(colour=NA))) + # removes black borders from legend
			coord_polar(theta='y') +
			theme(axis.ticks=element_blank(),
					axis.text.y=element_blank(),
					axis.text.x=element_blank(),
					axis.title=element_blank(),
					panel.grid=element_blank()) +
			scale_y_continuous(breaks=cumsum(datForFigure[[y]])-datForFigure[[y]]/2,
					labels=datForFigure[[fill]]) +
			theme(panel.background = element_rect(fill = "white"))+
			theme(legend.text=element_text(size=textSize),
					legend.title= element_text(size=textSize),
					strip.text.x = element_text(size = textSize))
	if (!is.na(main)) {
		p = p + ggtitle(main)
	}
	if (!is.na(facet)) {
		p<-p+facet_wrap(c(facet))
	}
	if (!is.na(colorNames)) {
		colors<-makeColors(length(unique(datForFigure[,fill])),colorNames)
		p<-p+scale_fill_manual(values=colors)
	}
	return(p)
}

ggpieToFile<-function(dat,fileName,fill="Category", maxCategory=5,textSize=9,transformTable=TRUE,...) {
	png(fileName,width=2000,height=2000,res=300)
	p<-ggpie(dat,fill=fill, maxCategory=maxCategory,textSize=textSize,transformTable=transformTable,...)
	print(p)
	dev.off()
}

ggpieGroupToFile<-function(dat,fileName,groupFileList="",outFileName="",
		maxCategory=5,textSize=9,transformTable=TRUE,fill="Category", y="Reads",facet="Sample",...) {
	if (groupFileList!="") {
		if (!transformTable) {
			dat<-acast(dat,as.formula(paste(fill,"~",facet)),value.var=y)
		}
		
		sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
		#keep the groups with samples in the count table
		sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(dat)),]
		
		datBySampleGroup<-mergeTableBySampleGroup(dat,sampleToGroup)
		if (outFileName!="") {
			write.csv(datBySampleGroup,outFileName)
		}
		ggpieToFile(datBySampleGroup,fileName=fileName,maxCategory=maxCategory,textSize=textSize,
				transformTable=TRUE,fill=fill,y=y,facet=facet,...)
	}
}

expandCountTableByName<-function(x,sep=";") {
	namesEach<-strsplit(row.names(x),sep)
	namesEachUnlist<-unlist(namesEach)
	namesEachLength<-sapply(namesEach,length)
	xEach<-x/namesEachLength
	result<-xEach[rep(seq.int(1,nrow(xEach)), namesEachLength),]
	row.names(result)<-namesEachUnlist
	return(result)
}

aggregateCountTable<-function(x,group,method=sum) {
	result<-aggregate(x, list(factor(group,levels=unique(group))), method)
	row.names(result)<-result[,1]
	result<-result[,-1]
	return(result)
}

countTableToSpecies<-function(dat,databaseLogFile="",outFileName="",shortName=T) {
	if (databaseLogFile!="") { #Species Count Table
		databaseLog<-read.delim(databaseLogFile,header=T,as.is=T)
		id2Species<-databaseLog$Species
		names(id2Species)<-databaseLog$Id
		
		mappingResultExpand<-expandCountTableByName(dat)
		speciesInMappingResult<-id2Species[row.names(mappingResultExpand)]
		
		mappingResult2Species<-aggregateCountTable(mappingResultExpand,speciesInMappingResult)
		
		#short name
		if (shortName) {
			row.names(mappingResult2Species)<-sapply(strsplit(row.names(mappingResult2Species),"_"),function(x) {
						if (length(x)<=3) {
							paste(x,collapse="_")
						} else if (grepl("^\\d+$",x[2])) {
							paste(x[1:3],collapse="_")
						} else {
							paste(x[1:2],collapse="_")
						}
					})
		}
		if (outFileName!="") {
			write.csv(mappingResult2Species,outFileName)
		}
	} else {
		mappingResult2Species<-dat
	}
	return(mappingResult2Species)
}

###############################################################################
# End funtions in count table barplot and pie chart
# Start defining parameters for functions
###############################################################################

if (!exists("maxCategory")) {
	maxCategory=5
}
if (!exists("textSize")) {
	textSize=9
}

###############################################################################
# End defining parameters for functions
###############################################################################
