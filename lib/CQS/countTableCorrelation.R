countTableFileList<-parSampleFile1
fixColorRange<-TRUE

library(heatmap3)
library(DESeq2)  
library(RColorBrewer)

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
	if (min(countNumCor,na.rm=T)<0) {
		if (fixColorRange) {
			col<-col_part(data_all=c(-1,1),data_part=countNumCor,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
		} else {
			col<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
		}
	} else {
		if (fixColorRange) {
			col<-col_part(data_all=c(0,1),data_part=countNumCor,col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))
		} else {
			col<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
		}
	}
	
	png(paste0(countTableFile,".Correlation.png"),width=2000,height=2000,res=300)
	heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,Rowv=NA,Colv=NA,col=col)
	dev.off()
	png(paste0(countTableFile,".Correlation.Cluster.png"),width=2000,height=2000,res=300)
	heatmap3(countNumCor,scale="none",balanceColor=T,margin=margin,col=col)
	dev.off()
}

