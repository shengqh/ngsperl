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

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

drawHCA<-function(prefix, rldselect, ispaired, designData, conditionColors, gnames){
  htfile<-paste0(prefix, "_DESeq2-vsd-heatmap.png")
  if(nrow(rldselect) > 2){
    png(filename=htfile, width=3000, height =3000, res=300)
    cexCol = max(1.0, 0.2 + 1/log10(ncol(rldselect)))
    if(ispaired){
      htColors<-rainbow(length(unique(designData$Paired)))
      gsColors<-as.matrix(data.frame(Group=conditionColors, Sample=htColors[designData$Paired]))
    }else{
      gsColors = conditionaColors;
    }
    heatmap3(rldselect, col = hmcols, ColSideColors = gsColors, margins=c(12,5), scale="r", dist=dist, labRow="",
    				 main=paste0("Hierarchical Cluster Using ", nrow(rldselect), " Genes"),  cexCol=cexCol, 
             legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
    dev.off()
  }
}

drawPCA<-function(prefix, rldmatrix, showLabelInPCA, designData){
  png(filename=paste0(prefix, "_DESeq2-vsd-pca.png"), width=3000, height=3000, res=300)
  pca<-prcomp(t(rldmatrix))
  supca<-summary(pca)$importance
  pcadata<-data.frame(pca$x)
  pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
  
  if(showLabelInPCA){
    g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=row.names(pcadata))) + 
      geom_text(vjust=-0.6, size=4) +
      geom_point(col=conditionColors, size=4) + 
      scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
      geom_hline(aes(0), size=.2) + 
      geom_vline(aes(0), size=.2) + 
      xlab(pcalabs[1]) + ylab(pcalabs[2])
  }else{
    g <- ggplot(pcadata, aes(x=PC1, y=PC2, color=designData$Condition)) + 
      geom_point(size=4) + 
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

countData<-data
index<-1
indecies<-c()
while(! is.numeric(countData[1,1])){
  countData<-countData[,c(2:ncol(countData))]
  indecies<-c(indecies, index)
  index<-index+1
}
countData[is.na(countData)] <- 0
countData<-round(countData)

comparisonNames=names(comparisons)
comparisonName=comparisonNames[1]

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
  
  if(ispaired){
    dir.create("spearman", showWarnings = FALSE)
    
    pairedSamples = unique(designData$Paired)
    
    spcorr<-unlist(lapply(c(1:length(pairedSamples)), function(x){
      samples<-designData$Sample[designData$Paired==pairedSamples[x]]
              cor(comparisonData[,samples[1]],comparisonData[,samples[2]],method="spearman")
            }))
            

    sptable<-data.frame(Name=pairedSamples, Spcorr=spcorr)
    write.csv(sptable, file=paste0(comparisonName, "_Spearman.csv"), row.names=FALSE)
    
    lapply(c(1:length(pairedSamples)), function(x){
      samples<-designData$Sample[designData$Paired==pairedSamples[x]]
      log2c1<-log2(comparisonData[,samples[1]]+1)
      log2c2<-log2(comparisonData[,samples[2]]+1)
      png(paste0("spearman/", comparisonName, "_Spearman_", pairedSamples[x], ".png"), width=2000, height=2000, res=300)
      plot(log2c1, log2c2, xlab=paste0(samples[1], " [log2(Count + 1)]"), ylab=paste0(samples[2], " [log2(Count + 1)]"))
      text(3,15,paste0("SpearmanCorr=", sprintf("%0.3f", spcorr[x])))
      dev.off()
    })
    
    pairedspearman[[comparisonName]]<-spcorr
  }
  
  notEmptyData<-apply(comparisonData, 1, max) > 0
  comparisonData<-comparisonData[notEmptyData,]
  
  if(ispaired){
    colnames(comparisonData)<-unlist(lapply(c(1:ncol(comparisonData)), function(i){paste0(designData$Paired[i], "_", colnames(comparisonData)[i])}))
  }
  rownames(designData)<-colnames(comparisonData)
  conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[designData$Condition]))
  
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
    inddata<-data[notEmptyData,indecies,drop=F]
    tbb<-cbind(inddata, comparisonData, res)
  }else{
    tbb<-cbind(comparisonData, res)
  }
  tbbselect<-tbb[select,,drop=F]
  
  tbb<-tbb[order(tbb$padj),,drop=F]
  write.csv(as.data.frame(tbb),paste0(comparisonName, "_DESeq2.csv"))
  
  tbbselect<-tbbselect[order(tbbselect$padj),,drop=F]
  write.csv(as.data.frame(tbbselect),paste0(comparisonName, "_DESeq2_sig.csv"))
  
  #some basic graph
  dds=DESeqDataSetFromMatrix(countData = comparisonData,
      colData = designData,
      design = ~1)
  
  colnames(dds)<-colnames(comparisonData)
  
  #draw density graph
  rldmatrix<-as.matrix(log2(counts(dds,normalized=FALSE) + 1))
  rsdata<-melt(rldmatrix)
  colnames(rsdata)<-c("Gene", "Sample", "log2Count")
  png(filename=paste0(comparisonName, "_DESeq2-log2-density.png"), width=4000, height=3000, res=300)
  g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + xlab("DESeq2 log2 transformed count")
  print(g)
  dev.off()
  
  #varianceStabilizingTransformation
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  assayvsd<-assay(vsd)
  write.csv(assayvsd, file=paste0(comparisonName, "_DESeq2-vsd.csv"))
  
  vsdiqr<-apply(assayvsd, 1, IQR)
  assayvsd<-assayvsd[order(vsdiqr, decreasing=T),]
  
  rldmatrix=as.matrix(assayvsd)

  siggenes<-rownames(rldmatrix) %in% rownames(tbbselect)
  
  #draw pca graph
  drawPCA(paste0(comparisonName,"_geneAll"), rldmatrix, showLabelInPCA, designData)
  if(showDEGeneCluster){
    drawPCA(paste0(comparisonName,"_geneNotDE"), rldmatrix[!siggenes,,drop=F], showLabelInPCA, designData)
  }
  
  #draw heatmap
  drawHCA(paste0(comparisonName,"_gene500"), rldmatrix[1:min(500, nrow(rldmatrix)),,drop=F], ispaired, designData, conditionColors, gnames)

  if(showDEGeneCluster){
	  drawHCA(paste0(comparisonName,"_geneDE"), , ispaired, designData, conditionColors, gnames)
	  nonDEmatrix<-rldmatrix[!siggenes,,drop=F]
	  drawHCA(paste0(comparisonName,"_geneAllNotDE"), nonDEmatrix, ispaired, designData, conditionColors, gnames)
	  drawHCA(paste0(comparisonName,"_gene500NotDE"), nonDEmatrix[1:min(500, nrow(rldmatrix)),,drop=F], ispaired, designData, conditionColors, gnames)
  }

  drawHCA(paste0(comparisonName,"_geneAll"), rldmatrix, ispaired, designData, conditionColors, gnames)
}

if(length(pairedspearman) > 0){
  #draw pca graph
  png(filename=paste0("spearman.png"), width=1000 * length(pairedspearman), height=2000, res=300)
  boxplot(pairedspearman)
  dev.off()
}
