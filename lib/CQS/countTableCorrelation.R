resultFile<-outFile
countTableFile<-parFile1

#resultFile<-"test"
#countTableFile<-"/scratch/cqs/zhaos/vickers/20160211_smallRNA_3018-KCV-52_53_54_human/bowtie1_genome_1mm_NTA_smallRNA_table/result/smallRNA_1mm_3018-KCV-52_53_54.count"

library(heatmap3)
library(DESeq2)  

if (grepl(".csv$",countTableFile)) {
	count<-read.csv(countTableFile,header=T,row.names=1)
} else {
	count<-read.delim(countTableFile,header=T,row.names=1)
}

colClass<-sapply(count, class)
countNum<-count[,which(colClass=="numeric" | colClass=="integer")]
countNum<-round(countNum,0)

dds=DESeqDataSetFromMatrix(countData = countNum, colData = as.data.frame(rep(1,ncol(countNum))),design = ~1)
temp<-DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
countNumVsd<-assay(temp)
colnames(countNumVsd)<-colnames(countNum)

countNumCor<-cor(countNumVsd,use="pa")

png(paste0(resultFile,".Correlation.png"),width=2000,height=2000,res=300)
heatmap3(countNumCor,scale="none",balanceColor=T,margin=c(8,8),Rowv=NA,Colv=NA)
dev.off()
png(paste0(resultFile,".Correlation.Cluster.png"),width=2000,height=2000,res=300)
heatmap3(countNumCor,scale="none",balanceColor=T,margin=c(8,8))
dev.off()
