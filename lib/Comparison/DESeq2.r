##predefined_condition_begin
#setwd("H:/shengquanhu/projects/smallRNA/20140926_bingshan_smallRNA_human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table_deseq2/result")  
#  
#data<-read.table("H:/shengquanhu/projects/smallRNA/20140926_bingshan_smallRNA_human/topN_bowtie1_genome_cutadapt_1mm_count_miRNA_table/result/miRNA_1mm_bingshan_smallRNA_human.count",row.names=1, header=T, check.names=F)
#
#comparisons=list(
#  "metastasis_vs_normal" = c("metastasis_vs_normal.design", "normal", "metastasis"),
#  "metastasis_vs_tumor" = c("metastasis_vs_tumor.design", "tumor", "metastasis"),
#  "tumor_vs_normal" = c("tumor_vs_normal.design", "normal", "tumor")
#)
##predefined_condition_end

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")

hmcols <- colorRampPalette(c("green", "black", "red"))(256)

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
  
#   if(ispaired){
#     pairedNames = unique(designData$Paired)
#     
#     spcorr<-unlist(lapply(c(1:length(g1)), function(x){
#               cor(c1[,x], c2[,x],method="spearman")
#             }))
#             
# 
#     sptable<-data.frame(Name=pairsamplenames, Spcorr=spcorr)
#     write.csv(sptable, file=paste0(pairname, "_Spearman.csv"), row.names=FALSE)
#     
#     lapply(c(1:length(g1)), function(x){
#       log2c1<-log2(c1[,x]+1)
#       log2c2<-log2(c2[,x]+1)
#       png(paste0(pairname, "_Spearman_", pairsamplenames[x], ".png"), width=2000, height=2000, res=300)
#       plot(log2c1, log2c2, xlab=paste0(g1[x], " [log2(Count + 1)]"), ylab=paste0(g2[x], " [log2(Count + 1)]"))
#       text(3,15,paste0("SpearmanCorr=", sprintf("%0.3f", cor(c1[,x], c2[,x],method="spearman")) ))
#       dev.off()
#     })
#     
#     pairedspearman[[pairname]]<-spcorr
#   }
  
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
  
  select<- (!is.na(res$padj)) & (res$padj<0.05) & ((res$log2FoldChange >= 1) | (res$log2FoldChange <= -1))
  
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
  
  #draw pca graph
  png(filename=paste0(comparisonName, "_DESeq2-vsd-pca.png"), width=3000, height=3000, res=300)
  pca<-prcomp(t(rldmatrix))
  supca<-summary(pca)$importance
  pcadata<-data.frame(pca$x)
  pcalabs=paste0(colnames(pcadata), "(", round(supca[2,] * 100), "%)")
  g <- ggplot(pcadata, aes(x=PC1, y=PC2, label=row.names(pcadata))) + 
      geom_text(vjust=-0.6, size=4) +
      geom_point(col=conditionColors, size=4) + 
      scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) +
      scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
      geom_hline(aes(0), size=.2) + 
      geom_vline(aes(0), size=.2) + 
      xlab(pcalabs[1]) + ylab(pcalabs[2])
  print(g)
  dev.off()
  
  #draw heatmap
  mincount<-min(500, nrow(rldmatrix))
  rldselect<-rldmatrix[1:mincount,,drop=F]
  htfile<-paste0(comparisonName, "_DESeq2-vsd-heatmap.png")
    if(nrow(rldselect) > 2){
      png(filename=htfile, width=3000, height =3000, res=300)
      if(ispaired){
        htColors<-rainbow(length(unique(designData$Paired)))
        gsColors<-as.matrix(data.frame(Group=conditionColors, Sample=htColors[designData$Paired]))
        heatmap3(rldselect, col = hmcols, ColSideColors = gsColors, margins=c(12,5), scale="r", dist=dist, labRow="", 
            legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
      }else{
        heatmap3(rldselect, col = hmcols, ColSideColors = pairColors, margins=c(12,5), scale="r", dist=dist, labRow="",
            legendfun=function() showLegend(legend=paste0("Group ", gnames),col=c("red","blue"),cex=1.0,x="center"))
      }
      dev.off()
    }
}
# 
# if(length(pairedspearman) > 0){
#   #draw pca graph
#   png(filename=paste0("spearman.png"), width=1000 * length(pairedspearman), height=2000, res=300)
#   boxplot(pairedspearman)
#   dev.off()
# }
