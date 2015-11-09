##predefined_condition_begin
#setwd("/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/star_deseq2/result")  
#  
#data<-read.table("/scratch/cqs/shengq1/rnaseq/20150226_bojana_FFPE_FF/hiseq/star_genetable/result/FFPE_FF_HiSeq_gene.count",row.names=1, header=T, check.names=F)
#
#showLabelInPCA<-1
#showDEGeneCluster<-1
#pvalue<-0.05
#foldChange<-2
#
#comparisons=list(
#  "HiSeq_FFPE2_VS_FF2_NoMismatch" = c("HiSeq_FFPE2_VS_FF2_NoMismatch.design", "HiSeq_FF_NoMismatch", "HiSeq_FFPE_NoMismatch")
#) 
#
##predefined_condition_end

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")
library("scales")

##Solving node stack overflow problem start###
#when there are too many genes, drawing dendrogram may failed due to node stack overflow,
#It could be solved by forcing stats:::plotNode to be run as interpreted code rather then byte-compiled code via a nasty hack.
#http://stackoverflow.com/questions/16559250/error-in-heatmap-2-gplots/25877485#25877485

# Convert a byte-compiled function to an interpreted-code function 
unByteCode <- function(fun)
{
  FUN <- eval(parse(text=deparse(fun)))
  environment(FUN) <- environment(fun)
  FUN
}

# Replace function definition inside of a locked environment **HACK** 
assignEdgewise <- function(name, env, value)
{
  unlockBinding(name, env=env)
  assign( name, envir=env, value=value)
  lockBinding(name, env=env)
  invisible(value)
}

# Replace byte-compiled function in a locked environment with an interpreted-code
# function
unByteCodeAssign <- function(fun)
{
  name <- gsub('^.*::+','', deparse(substitute(fun)))
  FUN <- unByteCode(fun)
  retval <- assignEdgewise(name=name,
                           env=environment(FUN),
                           value=FUN
  )
  invisible(retval)
}

# Use the above functions to convert stats:::plotNode to interpreted-code:
unByteCodeAssign(stats:::plotNode)

# Now raise the interpreted code recursion limit (you may need to adjust this,
#  decreasing if it uses to much memory, increasing if you get a recursion depth error ).
options(expressions=5e4)

##Solving node stack overflow problem end###

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

drawHCA<-function(prefix, rldselect, ispaired, designData, conditionColors, gnames){
  htfile<-paste0(prefix, "_DESeq2-vsd-heatmap.png")
  cat("saving HCA to ", htfile, "\n")
  genecount<-nrow(rldselect)
  if(genecount > 2){
    png(filename=htfile, width=3000, height =3000, res=300)
    cexCol = max(1.0, 0.2 + 1/log10(ncol(rldselect)))
    if(ispaired){
      htColors<-rainbow(length(unique(designData$Paired)))
      gsColors<-as.matrix(data.frame(Group=conditionColors, Sample=htColors[designData$Paired]))
    }else{
      gsColors = conditionColors;
    }
    heatmap3(rldselect, 
             col = hmcols, 
             ColSideColors = gsColors, 
             margins=c(12,5), 
             scale="r", 
             dist=dist, 
             labRow=NA,
    				 main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),  
             cexCol=cexCol, 
             legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
    dev.off()
  }
}

drawPCA<-function(prefix, rldmatrix, showLabelInPCA, designData, conditionColors){
  filename<-paste0(prefix, "_DESeq2-vsd-pca.png")
  genecount<-nrow(rldmatrix)
  if(genecount > 2){
	  cat("saving PCA to ", filename, "\n")
	  png(filename=filename, width=3000, height=3000, res=300)
	  pca<-prcomp(t(rldmatrix))
	  supca<-summary(pca)$importance
	  pcadata<-data.frame(pca$x)
	  pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
	  pcadata["sample"]<-row.names(pcadata)
	  
	  if(showLabelInPCA){
		  g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=sample)) + 
				  geom_text(vjust=-0.6, size=4) +
				  geom_point(col=conditionColors, size=4) + 
				  scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
				  scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
				  geom_hline(aes(0), size=.2) + 
				  geom_vline(aes(0), size=.2) + 
				  xlab(pcalabs[1]) + ylab(pcalabs[2])
	  }else{
		  g <- ggplot(pcadata, aes(x=PC1, y=PC2)) + 
				  geom_point(col=conditionColors, size=4) + 
				  labs(color = "Group") +
				  scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) + 
				  scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
				  geom_hline(aes(0), size=.2) + 
				  geom_vline(aes(0), size=.2) +
				  xlab(pcalabs[1]) + ylab(pcalabs[2]) + 
				  theme(legend.position="top")
	  }
	  
	  print(g)
	  dev.off()
  }
}

#for volcano plot
reverselog_trans <- function(base = exp(1)) {
	trans <- function(x) -log(x, base)
	inv <- function(x) base^(-x)
	trans_new(paste0("reverselog-", format(base)), trans, inv, 
			log_breaks(base = base), 
			domain = c(1e-100, Inf))
}

isDataNumeric = unlist(lapply(data[1,], function(x){is.numeric(x)}))
index = 1
while(!all(isDataNumeric[index:ncol(data)])){
  index = index + 1
}

if(index > 1){
  indecies<-c(1:(index-1))
}else{
  indecies<-c()
}
countData<-data[,c(index:ncol(data))]

countData[is.na(countData)] <- 0
countData<-round(countData)

comparisonNames=names(comparisons)
comparisonName=comparisonNames[1]

dir.create("details", showWarnings = FALSE)

pairedspearman<-list()
for(comparisonName in comparisonNames){
  str(comparisonName)
  designFile=comparisons[[comparisonName]][1]
  gnames=comparisons[[comparisonName]][2:3]
  designData<-read.table(designFile, sep="\t", header=T)
  designData$Condition<-factor(designData$Condition, levels=gnames)
  
  if(ncol(designData) == 3){
    ispaired<-TRUE
    cat("Paired data!\n")
  }else{
    ispaired<-FALSE
    cat("Not paired data!\n")
  }
  
  comparisonData<-countData[,colnames(countData) %in% as.character(designData$Sample),drop=F]
  if(ncol(comparisonData) != nrow(designData)){
    warning(paste0("Data not matched, there are ", nrow(designData), " samples in design file ", designFile, " but ", ncol(comparisonData), " samples in data "))
    next
  }
  comparisonData<-comparisonData[,as.character(designData$Sample)]
  
  prefix<-comparisonName
  curdata<-data
  if(minMedianInGroup > 0){
    conds<-unique(designData$Condition)
    data1<-comparisonData[, colnames(comparisonData) %in% designData$Sample[designData$Condition==conds[1]]]
    data2<-comparisonData[, colnames(comparisonData) %in% designData$Sample[designData$Condition==conds[2]]]
    med1<-apply(data1, 1, median) >= minMedianInGroup
    med2<-apply(data2, 1, median) >= minMedianInGroup
    med<-med1 | med2
    comparisonData<-comparisonData[med,]
    cat(nrow(comparisonData), " genes with minimum median count in group larger or equals than ", minMedianInGroup, "\n")
    prefix<-paste0(comparisonName, "_min", minMedianInGroup)
    curdata<-data[med,]
  }
  
  if(ispaired){
    pairedSamples = unique(designData$Paired)
    
    spcorr<-unlist(lapply(c(1:length(pairedSamples)), function(x){
      samples<-designData$Sample[designData$Paired==pairedSamples[x]]
              cor(comparisonData[,samples[1]],comparisonData[,samples[2]],method="spearman")
            }))
            

    sptable<-data.frame(Name=pairedSamples, Spcorr=spcorr)
    write.csv(sptable, file=paste0(prefix, "_Spearman.csv"), row.names=FALSE)
    
    lapply(c(1:length(pairedSamples)), function(x){
      samples<-designData$Sample[designData$Paired==pairedSamples[x]]
      log2c1<-log2(comparisonData[,samples[1]]+1)
      log2c2<-log2(comparisonData[,samples[2]]+1)
      png(paste0("details/", prefix, "_Spearman_", pairedSamples[x], ".png"), width=2000, height=2000, res=300)
      plot(log2c1, log2c2, xlab=paste0(samples[1], " [log2(Count + 1)]"), ylab=paste0(samples[2], " [log2(Count + 1)]"))
      text(3,15,paste0("SpearmanCorr=", sprintf("%0.3f", spcorr[x])))
      dev.off()
    })
    
    pairedspearman[[comparisonName]]<-spcorr
  }
  
  notEmptyData<-apply(comparisonData, 1, max) > 0
  comparisonData<-comparisonData[notEmptyData,]
  curdata<-curdata[notEmptyData,]
  
  if(ispaired){
    colnames(comparisonData)<-unlist(lapply(c(1:ncol(comparisonData)), function(i){paste0(designData$Paired[i], "_", colnames(comparisonData)[i])}))
  }
  rownames(designData)<-colnames(comparisonData)
  conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[designData$Condition]))
  
  #some basic graph
  dds=DESeqDataSetFromMatrix(countData = comparisonData,
      colData = designData,
      design = ~1)
  
  colnames(dds)<-colnames(comparisonData)
  
  #draw density graph
  rldmatrix<-as.matrix(log2(counts(dds,normalized=FALSE) + 1))
  rsdata<-melt(rldmatrix)
  colnames(rsdata)<-c("Gene", "Sample", "log2Count")
  png(filename=paste0(prefix, "_DESeq2-log2-density.png"), width=4000, height=3000, res=300)
  g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + xlab("DESeq2 log2 transformed count")
  print(g)
  dev.off()

  png(filename=paste0(prefix, "_DESeq2-log2-density-individual.png"), width=4000, height=3000, res=300)
	g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + facet_wrap(~Sample, scales = "free") + xlab("DESeq2 log2 transformed count")
  print(g)
  dev.off()
  
  #varianceStabilizingTransformation
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  assayvsd<-assay(vsd)
  write.csv(assayvsd, file=paste0(prefix, "_DESeq2-vsd.csv"))
  
  vsdiqr<-apply(assayvsd, 1, IQR)
  assayvsd<-assayvsd[order(vsdiqr, decreasing=T),]
  
  rldmatrix=as.matrix(assayvsd)

  #draw pca graph
  drawPCA(paste0(prefix,"_geneAll"), rldmatrix, showLabelInPCA, designData, conditionColors)
  
  #draw heatmap
  #drawHCA(paste0(prefix,"_gene500"), rldmatrix[1:min(500, nrow(rldmatrix)),,drop=F], ispaired, designData, conditionColors, gnames)
  drawHCA(paste0(prefix,"_geneAll"), rldmatrix, ispaired, designData, conditionColors, gnames)
  
  #different expression analysis
  if(ispaired){
    dds=DESeqDataSetFromMatrix(countData = comparisonData,
                               colData = designData,
                               design = ~ Paired + Condition)
  }else{
    dds=DESeqDataSetFromMatrix(countData = comparisonData,
                               colData = designData,
                               design = ~ Condition)
  }
  
  dds <- DESeq(dds)
  res<-results(dds,cooksCutoff=FALSE)
  
  cat("DESeq2 finished.\n")
  
  select<-(!is.na(res$padj)) & (res$padj<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))
  
  if(length(indecies) > 0){
    inddata<-curdata[,indecies,drop=F]
    tbb<-cbind(inddata, comparisonData, res)
  }else{
    tbb<-cbind(comparisonData, res)
  }
  tbb$FoldChange<-2^tbb$log2FoldChange
  tbbselect<-tbb[select,,drop=F]
  
  tbb<-tbb[order(tbb$padj),,drop=F]
  write.csv(as.data.frame(tbb),paste0(prefix, "_DESeq2.csv"))
  
  tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
  write.csv(as.data.frame(tbbselect),paste0(prefix, "_DESeq2_sig.csv"))
  
  if(showDEGeneCluster){
    siggenes<-rownames(rldmatrix) %in% rownames(tbbselect)

    nonDEmatrix<-rldmatrix[!siggenes,,drop=F]
    DEmatrix<-rldmatrix[siggenes,,drop=F]
    
    drawPCA(paste0(prefix,"_geneNotDE"), nonDEmatrix, showLabelInPCA, designData, conditionColors)
    drawHCA(paste0(prefix,"_geneNotDE"), nonDEmatrix, ispaired, designData, conditionColors, gnames)
    
	drawPCA(paste0(prefix,"_geneDE"),DEmatrix , showLabelInPCA, designData, conditionColors)
    drawHCA(paste0(prefix,"_geneDE"),DEmatrix , ispaired, designData, conditionColors, gnames)
    #drawHCA(paste0(prefix,"_gene500NotDE"), nonDEmatrix[1:min(500, nrow(nonDEmatrix)),,drop=F], ispaired, designData, conditionColors, gnames)
  }
  
  #Top 25 Significant genes barplot
  sigDiffNumber<-nrow(tbbselect)
  if (sigDiffNumber>0) {
	  if (sigDiffNumber>25) {
		  print(paste0("More than 25 genes were significant. Only the top 25 genes will be used in barplot"))
		  diffResultSig<-tbbselect[order(tbbselect$padj)[1:25],]
	  } else {
		  diffResultSig<-tbbselect
	  }
	  diffResultSig$Name<-sapply(strsplit(row.names(diffResultSig),";"),function(x) x[1])
	  diffResultSig$Name <- factor(diffResultSig$Name, levels=diffResultSig$Name[order(diffResultSig$log2FoldChange)])
	  diffResultSig<-as.data.frame(diffResultSig)
	  
	  png(filename=paste0(prefix, "_DESeq2_sig_barplot.png"), width=3000, height=3000, res=300)
#	  pdf(paste0(prefix,"_DESeq2_sig_barplot.pdf"))
	  p<-ggplot(diffResultSig,aes(x=Name,y=log2FoldChange,order=log2FoldChange))+geom_bar(stat="identity")+
			  coord_flip()+
#			geom_abline(slope=0,intercept=1,colour="red",linetype = 2)+
			  scale_y_continuous(name=bquote(log[2]~Fold~Change))+
			  theme(axis.text = element_text(colour = "black"))
	  print(p)
	  dev.off()
  } else {
	  print(paste0("No gene with adjusted p value less than ",pvalue," and fold change larger than ",foldChange))
  }
  
  #volcano plot
  changeColours<-c(grey="grey",green="green",red="red")
  diffResult<-as.data.frame(tbb)
  diffResult$log10BaseMean<-log10(diffResult$baseMean)
  diffResult$colour<-"grey"
  diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
  diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
  png(filename=paste0(prefix, "_DESeq2_volcanoPlot.png"), width=3000, height=3000, res=300)
#  pdf(paste0(prefix,"_DESeq2_volcanoPlot.pdf"))
  p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
		  geom_point(aes(size=log10BaseMean,colour=colour))+
		  scale_color_manual(values=changeColours,guide = FALSE)+
		  scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))+
		  scale_x_continuous(name=bquote(log[2]~Fold~Change))+
		  geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
		  geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
		  guides(size=guide_legend(title=bquote(log[10]~Base~Mean)))+
		  theme_bw()+
		  scale_size(range = c(2, 7))+
		  theme(axis.text = element_text(colour = "black"))
  print(p)
  dev.off()
}

if(length(pairedspearman) > 0){
  #draw pca graph
  filename<-ifelse(minMedianInGroup > 0, paste0("spearman_min", minMedianInGroup, ".png"), "spearman.png")
  png(filename=filename, width=1000 * length(pairedspearman), height=2000, res=300)
  boxplot(pairedspearman)
  dev.off()
}
