##predefined_condition_begin
#setwd("H:/shengquanhu/projects/chenxi/20131017_chenxi_rnaseq_smad4/deseq2/result")  
#data<-read.table("H:/shengquanhu/projects/chenxi/20131017_chenxi_rnaseq_smad4/genetable/result/smad4_gene.count",row.names=1, header=T, check.names=F)
#pairs=list(
#	"KO_vs_WT" = list("WT" = c("2288-RDB-81","2288-RDB-83","2288-RDB-85"), "KO" = c("2288-RDB-82","2288-RDB-84","2288-RDB-86")) 
#)
##predefined_condition_end

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

hasname <- (! is.numeric(data[1,1]))
if(hasname){
	countData<-data[,c(2:ncol(data))]
}else{
	countData<-data
}

pairnames=names(pairs)
pairname=pairnames[1]

for(pairname in pairnames){
	str(pairname)
	gs=pairs[[pairname]]
	gnames=names(gs)
	g1name=gnames[1]
	g2name=gnames[2]
	g1=gs[[g1name]]
	g2=gs[[g2name]]
	c1=countData[,colnames(countData) %in% g1,drop=F]
	c2=countData[,colnames(countData) %in% g2,drop=F]
	
	if(ncol(c1) != length(g1)){
		warning(paste0("There are only ", ncol(c1), " samples in group ", g1name, " but ", length(g1), " required!"))
		next
	}
	
	if(ncol(c2) != length(g2)){
		warning(paste0("There are only ", ncol(c2), " samples in group ", g2name, " but ", length(g2), " required!"))
		next
	}
	
	pairCountData=cbind(c1, c2)
	
	pairColData=data.frame(condition=factor(c(rep(g1name, ncol(c1)), rep(g2name, ncol(c2))), levels=gnames))
	rownames(pairColData)<-colnames(pairCountData)
	pairColors<-as.matrix(data.frame(Group=c(rep("red", ncol(c1)), rep("blue", ncol(c2)))))
	
	#different expression analysis
	dds=DESeqDataSetFromMatrix(countData = pairCountData,
			colData = pairColData,
			design = ~ condition)
	
	dds <- DESeq(dds)
	res<-results(dds,cooksCutoff=FALSE)
	
	select<- (!is.na(res$padj)) & (res$padj<0.05) & ((res$log2FoldChange >= 1) | (res$log2FoldChange <= -1))
	
	if(hasname){
		tbb<-cbind(data[,1,drop=F], pairCountData, res)
	}else{
		tbb<-cbind(pairCountData, res)
	}
	tbbselect<-tbb[select,,drop=F]
	
	tbb<-tbb[order(tbb$padj),,drop=F]
	write.csv(as.data.frame(tbb),paste0(pairname, "_DESeq2.csv"))
	
	tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
	write.csv(as.data.frame(tbbselect),paste0(pairname, "_DESeq2_sig.csv"))
	
	rld<-rlogTransformation(dds, blind=TRUE)
	rldmatrix=as.matrix(assay(rld))
	
	#draw density graph
	rsdata<-melt(rldmatrix)
	colnames(rsdata)<-c("Gene", "Sample", "logCount")
	png(filename=paste0(pairname, ".logdensity.png"), width=4000, height=3000, res=300)
	g<-ggplot(rsdata) + geom_density(aes(x=logCount, colour=Sample)) + xlab("DESeq2 log transformed count")
	print(g)
	dev.off()
	
	#draw pca graph
	png(filename=paste0(pairname, ".pca.png"), width=3000, height=3000, res=300)
	pca<-prcomp(t(rldmatrix))
	supca<-summary(pca)$importance
	data<-data.frame(pca$x)
	pcalabs=paste0(colnames(data), "(", round(supca[2,] * 100), "%)")
	g <- ggplot(data, aes(x=PC1, y=PC2, label=row.names(data))) + 
			geom_text(vjust=-0.6, size=4) +
			geom_point(col=pairColors, size=4) + 
			scale_x_continuous(limits=c(min(data$PC1) * 1.2,max(data$PC1) * 1.2)) +
			scale_y_continuous(limits=c(min(data$PC2) * 1.2,max(data$PC2) * 1.2)) + 
			geom_hline(aes(0), size=.2) + 
			geom_vline(aes(0), size=.2) + 
			xlab(pcalabs[1]) + ylab(pcalabs[2])
	print(g)
	dev.off()
	
	#draw heatmap
	rldselect<-rldmatrix[select,,drop=F]
	if(nrow(rldselect) > 2){
		png(filename=paste0(pairname, ".heatmap.png"), width=3000, height =3000, res=300)
		heatmap3(rldselect, col = hmcols, ColSideColors = pairColors, margins=c(12,5), scale="r", dist=dist, labRow="",
				legendfun=function() showLegend(legend=paste0("Group ", gnames),col=c("red","blue"),cex=1.5,x="center"))
		dev.off()
	}
}