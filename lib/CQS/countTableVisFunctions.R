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

tablePie<-function(x,maxCategory=10,main="",addPercent=F) {
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
	if (n<8) {
		colors<-brewer.pal(n, colorNames)
		return(colors)
	} else {
		
	}
}

tableMaxCategory<-function(dat,maxCategory=NA) {
	if (!is.na(maxCategory) & nrow(dat)>maxCategory) {
		temp<-apply(dat,2,function(x) rev(order(x))[1:maxCategory])
		categoryKeptInd<-sort(unique(as.vector(temp)))
		datForFigure<-dat[categoryKeptInd,]
		if (length(categoryKeptInd)<nrow(dat)) {
			datForFigure<-rbind(datForFigure,Other=colSums(dat[-categoryKeptInd,]))
		}
	} else {
		datForFigure<-dat
	}
	return(datForFigure)
}

tableBarplot<-function(dat,maxCategory=5,x="Sample", y="Reads",fill="Category",groupName=fill,transformTable=T,textSize=20) {
	if (transformTable) {
		datForFigure<-tableMaxCategory(dat,maxCategory=maxCategory)
		
		datForFigure$Groups<-row.names(dat)
		datForFigure<-melt(xForFigure)
		colnames(xForFigure)<-c(fill,x,y)
	} else {
		datForFigure<-dat
	}
	p<-ggplot(xForFigure,aes_string(x=x,y=y,fill=fill))+
			geom_bar(stat="identity")+
			guides(fill= guide_legend(title = groupName))+
			theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
			theme(legend.position = "top")+
			theme(axis.text = element_text(size=textSize),legend.text=element_text(size=textSize),
					axis.title = element_text(size=textSize),legend.title= element_text(size=textSize))+
			guides(fill = guide_legend(keywidth = 2, keyheight = 2))
	return(p)
}

#changed from function in http://mathematicalcoffee.blogspot.com/2014/06/ggpie-pie-graphs-in-ggplot2.html
ggpie <- function (dat, fill="Species", y="Reads",facet="Sample", maxCategory=NA,main=NA, percent=T,textSize=15) {
	datForFigure<-tableMaxCategory(dat,maxCategory=maxCategory)

	if (percent) {
		datForFigure<-t(t(datForFigure)/colSums(datForFigure))
	}
	if (row.names(datForFigure)[nrow(datForFigure)]=="Other") {
		categoryOrderedNames<-c(row.names(datForFigure)[-nrow(datForFigure)][rev(order(rowSums(datForFigure[-nrow(datForFigure),])))],"Other")
	} else {
		categoryOrderedNames<-row.names(datForFigure)[rev(order(rowSums(datForFigure)))]
	}
	
	datForFigure<-melt(as.matrix(datForFigure),as.is=T)
	colnames(datForFigure)<-c(fill,facet,y)
	datForFigure[,1]<-factor(datForFigure[,1],levels=categoryOrderedNames)
	datForFigure<-orderDataByNames(datForFigure,datForFigure[,1],categoryOrderedNames)
	
	p = ggplot(datForFigure, aes_string(x=factor(1), y=y, fill=fill)) +
			geom_bar(width=1, stat='identity', color='black') +
			guides(fill=guide_legend(keywidth = 1.5, keyheight = 1.5,override.aes=list(colour=NA))) + # removes black borders from legend
			coord_polar(theta='y') +
			theme(axis.ticks=element_blank(),
					axis.text.y=element_blank(),
					axis.text.x=element_blank(),
					axis.title=element_blank(),
					panel.grid=element_blank()) +
			scale_y_continuous(breaks=cumsum(datForFigure[[y]]) - datForFigure[[y]] / 2, labels=datForFigure[[fill]]) +
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
	return(p)
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



