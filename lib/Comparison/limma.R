
library(limma)



makeDesignTable <- function(rawDataTable, factorTable,factorInDesign=NULL) {
  if (is.null(factorInDesign)) {
    factorInDesign=as.character(unique(factorTable[, 2]))
  }
  designTable <- matrix(0, nrow = ncol(rawDataTable), ncol = length(factorInDesign))
  row.names(designTable) <- colnames(rawDataTable)
  colnames(designTable) <- factorInDesign
  
  designTable[as.matrix(factorTable)] <- 1
  designTable
  designTable <- cbind(Intercept = 1, designTable)
  return(designTable)
}



performLimma <- function(rawDataTable, designTable, contrastsText = NULL) {
  rawDataTable=rawDataTable[,row.names(designTable)]
  if (is.null(contrastsText)) {
    contrastsText <- c(colnames(designTable)[-1])
  }
  
  fit <- lmFit(rawDataTable, designTable)
  contrast.matrix <- makeContrasts(contrasts = contrastsText, levels = designTable)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  return(fit2)
}



performNormlization = TRUE
pvalue=0.05
useRawPvalue=1
foldChange=0
notNaPropotion=0.5

log2FcCol="logFC"
statCol="t"
pvalueCol="P.Value"
adjPvalueCol="adj.P.Val"


########################################################
# Read and format data
########################################################

rawDataFile <- parFile1
rawData <- read.csv(rawDataFile, header = T, row.names = 1, as.is = T)

colClass <- sapply(rawData, class)
countNotNumIndex <- which((colClass != "numeric" & colClass != "integer") | grepl("Gene_Id", colnames(rawData)))
if (length(countNotNumIndex) == 0) {
  rawDataTable <- rawData
  rawDataAnnotation <- NULL
} else {
  rawDataTable <- rawData[, c(max(countNotNumIndex) + 1):ncol(rawData), drop = FALSE]
  rawDataAnnotation <- rawData[, c(1:max(countNotNumIndex)), drop = FALSE]
}


factorTable=read.delim(parSampleFile1,header=F,as.is=T)
factorInComparisons=read.delim(parSampleFile2,header=F,as.is=T)

################################################
# Normlization and boxplot/density plot before and after
################################################
png(paste0(outFile,".RawBoxplot.png"))
boxplot(rawDataTable, las = 2)
dev.off()

dataForPlot=melt(rawDataTable,variable.name="Sample")
p=ggplot(dataForPlot,aes(x=value,colour=Sample))+geom_density()
png(paste0(outFile,".RawDensity.png"))
plot(p)
dev.off()

if (performNormlization) {
  qua.norm <- function(pep) {
    require(affy)
    require(preprocessCore)
    pep.qua <- normalize.quantiles(as.matrix(pep))
    row.names(pep.qua) <- row.names(pep)
    colnames(pep.qua) <- colnames(pep)
    return(as.data.frame(pep.qua))
  }
  rawDataTable <- qua.norm(rawDataTable)
  png(paste0(outFile,".NormlizationBoxplot.png"))
  boxplot(rawDataTable, las = 2)
  dev.off()
  
  dataForPlot=melt(rawDataTable,variable.name="Sample")
  p=ggplot(dataForPlot,aes(x=value,colour=Sample))+geom_density()
  png(paste0(outFile,".NormlizationDensity.png"))
  plot(p)
  dev.off()
}




##################################
# heatmap and PCA
##################################
library(ComplexHeatmap)

designTableForHeatmap <- makeDesignTable(rawDataTable, factorTable)
designTableForHeatmap <- as.data.frame(designTableForHeatmap[,-1,drop=FALSE])
for (i in 1:ncol(designTableForHeatmap)) {
  designTableForHeatmap[[i]] <- as.character(designTableForHeatmap[[i]])
}
ha <- HeatmapAnnotation(df = designTableForHeatmap)
naCount=apply(rawDataTable,1,function(x) length(which(is.na(x))))
rawDataTableForHeatmap=rawDataTable[which(naCount/ncol(rawDataTable)<=1-notNaPropotion),]
ht <- Heatmap(as.matrix(rawDataTableForHeatmap), cluster_rows = FALSE, top_annotation = ha, row_names_gp = gpar(fontsize = 6))

png(paste0(outFile,".Heatmap.png"))
print(ht)
dev.off()

##################################
# For each comprasion, perform limma
##################################
naCount=apply(rawDataTable,1,function(x) length(which(is.na(x))))
rawDataTableForDiffDetection=rawDataTable[which(naCount/ncol(rawDataTable)<=1-notNaPropotion),]

DiffResultsOutTableAll=list() #record all results when needed
DiffResultsOutDiffGenesAll=list() #record all diff genes when needed
for (ComparisonOne in unique(factorInComparisons[,3])) {
#  prefix=paste0(outFile,"_",ComparisonOne)
  prefix=paste0(ComparisonOne)
  
  factorInDesign=unique(factorInComparisons[which(factorInComparisons[,3]==ComparisonOne),1])
  designTable <- makeDesignTable(rawDataTableForDiffDetection, factorTable,factorInDesign=factorInDesign)
  
  fit3 <- performLimma(rawDataTableForDiffDetection, designTable)
  DiffResults <- topTable(fit3, coef = colnames(designTable)[2],number=99999)
  
  ## Output differential genes table
  if (!is.null(rawDataAnnotation)) {
    tbb=cbind(rawDataAnnotation[row.names(DiffResults),],DiffResults)
  } else {
    tbb=DiffResults
  }
  if (useRawPvalue==1) {
    selectDiff<-(!is.na(tbb[,pvalueCol])) & (tbb[,pvalueCol]<pvalue) & ((tbb[,log2FcCol] >= log2(foldChange)) | (tbb[,log2FcCol] <= -log2(foldChange)))
  } else {
    selectDiff<-(!is.na(tbb[,adjPvalueCol])) & (tbb[,adjPvalueCol]<pvalue) & ((tbb[,log2FcCol] >= log2(foldChange)) | (tbb[,log2FcCol] <= -log2(foldChange)))
  }
  write.csv(as.data.frame(tbb),paste0(prefix, "_Diff.csv"))
  write.csv(as.data.frame(tbb[selectDiff,]),paste0(prefix, "_Diff_sig.csv"))
  
  ##################################
  # Output formated files for WebGesnet, GSEA, KEGGprofile
  ##################################
  geneNameField = NULL
  lowColNames = tolower(colnames(tbb))
  for(name in c("Feature_gene_name", "Gene.Symbol", "Gene_Symbol", "Gene Symbol", "GeneSymbol")){
    lowName = tolower(name)
    if(lowName %in% lowColNames){
      geneNameField=colnames(tbb)[match(lowName, lowColNames)]
      break
    }
  }
  
  if(!is.null(geneNameField)){
    write.table(na.omit(tbb[,c(geneNameField, statCol),drop=F]),paste0(prefix, "_GSEA.rnk"),row.names=F,col.names=F,sep="\t", quote=F)
    write.table(tbb[selectDiff,c(geneNameField),drop=F], paste0(prefix, "_sig_genename.txt"),row.names=F,col.names=F,sep="\t", quote=F)
  }else{
    warning("Can't identify Gene Name column, use row names as gene name")
    write.table(na.omit(tbb[,c(statCol),drop=F]),paste0(prefix, "_GSEA.rnk"),row.names=T,col.names=F,sep="\t", quote=F)
    write.table(data.frame(name=rownames(tbb[selectDiff,])), paste0(prefix, "_sig_genename.txt"),row.names=F,col.names=F,sep="\t", quote=F)
  }
  
  ##################################
  #save result for all
  ##################################
  DiffResultsOutTableAll[[ComparisonOne]]=tbb
  if (length(which(selectDiff))) {
    if (!is.null(geneNameField)) {
      DiffResultsOutDiffGenesAll[[ComparisonOne]]=na.omit(tbb[selectDiff,geneNameField])
    } else {
      DiffResultsOutDiffGenesAll[[ComparisonOne]]=na.omit(row.names(tbb)[selectDiff])
    }
  }
}


########################################################
# Export figures based on all results
########################################################
#vocalo plot
library(data.table)
dataForPlot=rbindlist(DiffResultsOutTableAll,idcol="Comparison")
if (useRawPvalue==1) {
  p=ggplot(dataForPlot,aes(x=logFC,y=-log10(P.Value)))
} else {
  p=ggplot(dataForPlot,aes(x=logFC,y=-log10(adj.P.Val)))
}
p=p+geom_point()+facet_wrap(~Comparison,scales="free")
png(paste0(outFile,".VolcanoPlot.png"))
print(p)
dev.off()

#venn plot of differential genes
#library("VennDiagram")
if (length(DiffResultsOutDiffGenesAll)>=1 & length(DiffResultsOutDiffGenesAll)<=5) {
  png(paste0(outFile,".Venn.png"))
  venn.diagram1(DiffResultsOutDiffGenesAll)
  dev.off()
}



save.image(paste0(outFile,".limma.RData"))














