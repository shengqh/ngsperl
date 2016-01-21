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

library("edgeR")

isDataNumeric = unlist(lapply(data[1,], function(x){is.numeric(x)}))
index = 1
while(!all(isDataNumeric[index:ncol(data)])){
  index = index + 1
}

indecies<-c(1:(index-1))
countData<-data[,c(index:ncol(data))]

countData[is.na(countData)] <- 0
countData<-round(countData)

comparisonNames=names(comparisons)
comparisonName=comparisonNames[1]

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
  
  notEmptyData<-apply(comparisonData, 1, max) > 0
  comparisonData<-comparisonData[notEmptyData,]
  curdata<-curdata[notEmptyData,]
  
  if(ispaired){
    colnames(comparisonData)<-unlist(lapply(c(1:ncol(comparisonData)), function(i){paste0(designData$Paired[i], "_", colnames(comparisonData)[i])}))
  }
  rownames(designData)<-colnames(comparisonData)
  conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[designData$Condition]))
  
  #different expression analysis
	condition<-designData$Condition
  if(ispaired){
  	subject<-designData$Paired
  	design <- model.matrix(~subject+condition)
  }else{
  	design <- model.matrix(~condition)
  }
  y<-DGEList(comparisonData)
  y <- estimateGLMCommonDisp(y,design)
  y <- estimateGLMTrendedDisp(y,design)
  y <- estimateGLMTagwiseDisp(y,design)
  fit <- glmFit(y,design)
  lrt <- glmLRT(fit,coef=2)
  res<-topTags(lrt)
    
  cat("EdgeR finished.\n")
  
  select<-(!is.na(res$padj)) & (res$padj<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))
  
  if(length(indecies) > 0){
    inddata<-curdata[,indecies,drop=F]
    tbb<-cbind(inddata, comparisonData, res)
  }else{
    tbb<-cbind(comparisonData, res)
  }
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
    
    drawHCA(paste0(prefix,"_geneDE"),DEmatrix , ispaired, designData, conditionColors, gnames)
    #drawHCA(paste0(prefix,"_gene500NotDE"), nonDEmatrix[1:min(500, nrow(nonDEmatrix)),,drop=F], ispaired, designData, conditionColors, gnames)
  }
}

if(length(pairedspearman) > 0){
  #draw pca graph
  filename<-ifelse(minMedianInGroup > 0, paste0("spearman_min", minMedianInGroup, ".png"), "spearman.png")
  png(filename=filename, width=1000 * length(pairedspearman), height=2000, res=300)
  boxplot(pairedspearman)
  dev.off()
}
