countTableFileList<-parSampleFile1

library(heatmap3)
library(DESeq2)  

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
	
	countNumCor<-cor(countNumVsd,use="pa")
	if (min(countNumCor,na.rm=T)<=0) {
		col<-colorRampPalette(c("navy", "white","firebrick3"))(1024)
	} else {
		col<-colorRampPalette(c("white","firebrick3"))(1024)
	}
	
	png(paste0(countTableFile,".Correlation.png"),width=2000,height=2000,res=300)
	heatmap3(countNumCor,scale="none",balanceColor=T,margin=c(8,8),Rowv=NA,Colv=NA,col=col)
	dev.off()
	png(paste0(countTableFile,".Correlation.Cluster.png"),width=2000,height=2000,res=300)
	heatmap3(countNumCor,scale="none",balanceColor=T,margin=c(8,8),col=col)
	dev.off()
}

