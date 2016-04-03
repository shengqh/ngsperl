countTableFileList<-parSampleFile1
groupFileList<-parSampleFile2
fixColorRange<-TRUE

library(heatmap3)
library(DESeq2)  
library(RColorBrewer)

#Merge Sample By Group
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

#extract part of color from a color range
col_part<-function(data_all,data_part,col) {
	min_all<-min(data_all,na.rm=T)
	max_all<-max(data_all,na.rm=T)
	min_part<-min(data_part,na.rm=T)
	max_part<-max(data_part,na.rm=T)
	cut_off_low<-round(quantile(1:length(col),(min_part-min_all)/(max_all-min_all)))
	cut_off_high<-round(quantile(1:length(col),(max_part-min_all)/(max_all-min_all)))
	col=col[cut_off_low:cut_off_high]
	return(col)
}

countTableFileAll<-read.delim(countTableFileList,header=F,as.is=T)
for (i in 1:nrow(countTableFileAll)) {
	countTableFile<-countTableFileAll[i,1]
	
	if (grepl(".csv$",countTableFile)) {
		count<-read.csv(countTableFile,header=T,row.names=1,as.is=T)
	} else {
		count<-read.delim(countTableFile,header=T,row.names=1,as.is=T)
	}
	
	colClass<-sapply(count, class)
	countNum<-count[,which(colClass=="numeric" | colClass=="integer")]
	countNum<-round(countNum,0)
	#remove genes with total reads 0
	countNum<-countNum[which(rowSums(countNum,na.rm=T)>0),]
	
	dds=DESeqDataSetFromMatrix(countData = countNum, colData = as.data.frame(rep(1,ncol(countNum))),design = ~1)
	temp<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
	countNumVsd<-assay(temp)
	colnames(countNumVsd)<-colnames(countNum)
	
	#heatmap
	margin=c(min(10,max(nchar(colnames(countNumVsd)))/2),min(10,max(nchar(row.names(countNumVsd)))/2))
	png(paste0(countTableFile,".heatmap.png"),width=2000,height=2000,res=300)
	heatmap3(countNumVsd,dist=dist,margin=margin,balanceColor=TRUE,,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
	dev.off()
	
	#correlation distribution
	countNumCor<-cor(countNumVsd,use="pa",method="sp")
	margin=c(min(10,max(nchar(colnames(countNumCor)))/2),min(10,max(nchar(row.names(countNumCor)))/2))
	
	colAll<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
	if (min(countNumCor,na.rm=T)<0) {
		colAllLabel<-c(-1,0,1)
		if (fixColorRange) {
			col<-col_part(data_all=c(-1,1),data_part=countNumCor,col=colAll)
		} else {
			col<-colAll
		}
	} else {
		colAllLabel<-c(0,0.5,1)
		if (fixColorRange) {
			col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colAll)
		} else {
			col<-colAll
		}
	}
	
	legendfun<-function(x) {
		par(mar = c(5, 1, 1, 1));
		image(x=1:length(colAll),y=1,z=matrix(1:length(colAll),ncol=1),xlab="",xaxt="n",yaxt="n",col=colAll);
		axis(1,at=c(1,length(colAll)/2,length(colAll)),labels=colAllLabel)
	}
	
	png(paste0(countTableFile,".Correlation.png"),width=2000,height=2000,res=300)
	heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
	dev.off()
	png(paste0(countTableFile,".Correlation.Cluster.png"),width=2000,height=2000,res=300)
	heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col,legendfun=legendfun)
	dev.off()
	
	
	if (groupFileList!="") {
		sampleToGroup<-read.delim(groupFileList,as.is=T,header=F)
		#keep the groups with samples in the count table
		sampleToGroup<-sampleToGroup[which(sampleToGroup[,1] %in% colnames(countNumVsd)),]
		countNumVsdGroup<-mergeTableBySampleGroup(countNumVsd,sampleToGroup)
		
		#heatmap
		margin=c(min(10,max(nchar(colnames(countNumVsdGroup)))/2),min(10,max(nchar(row.names(countNumVsdGroup)))/2))
		png(paste0(countTableFile,".Group.heatmap.png"),width=2000,height=2000,res=300)
		heatmap3(countNumVsdGroup,dist=dist,margin=margin,balanceColor=TRUE,,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
		dev.off()
		
		#correlation distribution
		countNumCor<-cor(countNumVsdGroup,use="pa",method="sp")
		margin=c(min(10,max(nchar(colnames(countNumCor)))/2),min(10,max(nchar(row.names(countNumCor)))/2))
		
		colAll<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
		if (min(countNumCor,na.rm=T)<0) {
			colAllLabel<-c(-1,0,1)
			if (fixColorRange) {
				col<-col_part(data_all=c(-1,1),data_part=countNumCor,col=colAll)
			} else {
				col<-colAll
			}
		} else {
			colAllLabel<-c(0,0.5,1)
			if (fixColorRange) {
				col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colAll)
			} else {
				col<-colAll
			}
		}
		
		legendfun<-function(x) {
			par(mar = c(5, 1, 1, 1));
			image(x=1:length(colAll),y=1,z=matrix(1:length(colAll),ncol=1),xlab="",xaxt="n",yaxt="n",col=colAll);
			axis(1,at=c(1,length(colAll)/2,length(colAll)),labels=colAllLabel)
		}
		
		
		png(paste0(countTableFile,".Group.Correlation.png"),width=2000,height=2000,res=300)
		heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col,legendfun=legendfun)
		dev.off()
		png(paste0(countTableFile,".Group.Correlation.Cluster.png"),width=2000,height=2000,res=300)
		heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col,legendfun=legendfun)
		dev.off()
		
	}
	
}



