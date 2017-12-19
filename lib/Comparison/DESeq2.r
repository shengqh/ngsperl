##predefined_condition_begin

rootdir<-"/scratch/cqs/shengq1/vickers/20170222_smallRNA_3018_61_human_v3/host_genome/deseq2_miRNA/result"
inputfile<-"3018_61.define" 

showLabelInPCA<-1
showDEGeneCluster<-1
pvalue<-0.05
foldChange<-1.5
minMedianInGroup<-5
addCountOne<-0
usePearsonInHCA<-0
top25only<-0
detectedInBothGroup<-1
performWilcox<-0
useRawPvalue<-1
textSize<-9
transformTable<-0
exportSignificantGeneName<-1

libraryFile<-"/scratch/cqs/shengq1/vickers/20170222_smallRNA_3018_61_human_v3/host_genome/bowtie1_genome_1mm_NTA_smallRNA_category/result/3018_61.Category.Table.csv"
libraryKey<-"TotalReads"

##predefined_condition_end

options(bitmapType='cairo')

suffix<-"";
if(top25only){
  suffix=paste0(suffix,"_top25")
}

if(detectedInBothGroup){
  suffix=paste0(suffix, "_detectedInBothGroup")
}

if(minMedianInGroup > 0){
  suffix=paste0(suffix, "_min", minMedianInGroup)
}

zeroCount=0
if(addCountOne){
  zeroCount=1
  minMedianInGroup=minMedianInGroup+1
}

if(!exists("usePearsonInHCA")){
  usePearsonInHCA=0
}

if(!exists("exportSignificantGeneName")){
  exportSignificantGeneName<-1
}

if(exists("libraryFile")){
  librarySize<-read.csv(libraryFile, row.names=1,check.names=FALSE)
  librarySize<-unlist(librarySize[libraryKey,,drop=T])
  cat("Using ", libraryKey, " in " , libraryFile , " as library size. \n")
}

library("DESeq2")
library("heatmap3")
library("lattice")
library("reshape")
library("ggplot2")
library("grid")
library("scales")
library("reshape2")
library("VennDiagram")
library("RColorBrewer")
library("preprocessCore")

setwd(rootdir)  
comparisons_data<-read.table(inputfile, header=T, check.names=F , sep="\t", stringsAsFactors = F)

##Solving node stack overflow problem start###
#when there are too many genes, drawing dendrogram may failed due to node stack overflow,
#It could be solved by forcing stats:::plotNode to be run as interpreted code rather then byte-compiled code via a nasty hack.
#http://stackoverflow.com/questions/16559250/error-in-heatmap-2-gplots/25877485#25877485

#align two count table
align<-function(data1,data2,by=0,suffixes=c(deparse(substitute(data1)),deparse(substitute(data2))),sort=T) {
  if (is.null(data1)) {
    return(data2)
  } else if (is.null(data2)) {
    return(data1)
  }
  data<-merge(data1,data2,by=by,all=T,suffixes=suffixes,sort=sort)
  row.names(data)<-data[,1]
  data<-data[,-1]
  return (data)
}

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
  showRowDendro = genecount <= 50
  if(genecount > 2){
    png(filename=htfile, width=3000, height=3000, res=300)
    cexCol = max(1.0, 0.2 + 1/log10(ncol(rldselect)))
    if(ispaired){
      htColors<-rainbow(length(unique(designData$Paired)))
      gsColors<-as.matrix(data.frame(Group=conditionColors, Sample=htColors[designData$Paired]))
    }else{
      gsColors = conditionColors;
    }
    if(usePearsonInHCA){
      heatmap3(rldselect, 
               col = hmcols, 
               ColSideColors = gsColors, 
               margins=c(12,5), 
               scale="r", 
               labRow=NA,
               showRowDendro=showRowDendro,
               main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),  
               cexCol=cexCol, 
               useRaster=FALSE,
               legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
    }else{
      heatmap3(rldselect, 
               col = hmcols, 
               ColSideColors = gsColors, 
               margins=c(12,5), 
               scale="r", 
               distfun=dist, 
               labRow=NA,
               showRowDendro=showRowDendro,
               main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),  
               cexCol=cexCol, 
               useRaster=FALSE,
               legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
    }
    dev.off()
  }
}

drawPCA<-function(prefix, rldmatrix, showLabelInPCA, designData, conditionColors){
  #filename<-paste0(prefix, "_DESeq2-vsd-pca.png")
  filename<-paste0(prefix, "_DESeq2-vsd-pca.pdf")
  genecount<-nrow(rldmatrix)
  if(genecount > 2){
    cat("saving PCA to ", filename, "\n")
    #png(filename=filename, width=3000, height=3000, res=300) # 10 X 10 inches
    pdf(filename, width=10, height=10)
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
        geom_hline(aes(yintercept=0), size=.2) + 
        geom_vline(aes(xintercept=0), size=.2) + 
        xlab(pcalabs[1]) + ylab(pcalabs[2])
    }else{
      g <- ggplot(pcadata, aes(x=PC1, y=PC2)) + 
        geom_point(col=conditionColors, size=4) + 
        labs(color = "Group") +
        scale_x_continuous(limits=c(min(pcadata$PC1) * 1.2,max(pcadata$PC1) * 1.2)) + 
        scale_y_continuous(limits=c(min(pcadata$PC2) * 1.2,max(pcadata$PC2) * 1.2)) + 
        geom_hline(aes(yintercept=0), size=.2) + 
        geom_vline(aes(xintercept=0), size=.2) +
        xlab(pcalabs[1]) + ylab(pcalabs[2]) + 
        theme(legend.position="top")
    }
    
    print(g)
    dev.off()
  }
}

myEstimateSizeFactors<-function(dds){
  if(exists("librarySize")){
    curLibrarySize<-librarySize[colnames(dds)]
    curSizeFactor<-curLibrarySize/median(curLibrarySize)
    sizeFactors(dds)<-curSizeFactor
  }else{
    sfres<-try(dds<-estimateSizeFactors(dds))
    if (class(sfres) == "try-error") {
      library(edgeR)
      y<-DGEList(counts=counts(dds))
      y<-calcNormFactors(y, methold="TMM")
      sizeFactors(dds)<-y$samples$norm.factors
    }
  }
  return(dds)
}

#for volcano plot
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

countfiles<-unlist(unique(comparisons_data$CountFile))
allComparisons<-unlist(unique(comparisons_data$ComparisonName))
if(length(allComparisons) != nrow(comparisons_data)){
  error(paste("Comparison names cannot be repeated ", comparisons_data$ComparisonName, sep=": "))
}
allTitles<-comparisons_data$ComparisonTitle
names(allTitles)<-comparisons_data$ComparisonName

dataAllOut<-NULL
resultAllOut<-NULL

allSigNameList<-list()
allSigDirectionList<-list()
sigTableAll<-NULL
sigTableAllGene<-NULL
sigTableAllVar<-c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","FoldChange")

countfile_index = 1
titles<-NULL
validComparisons<-c()
for(countfile_index in c(1:length(countfiles))){
  countfile = countfiles[countfile_index]
  comparisons = comparisons_data[comparisons_data$CountFile == countfile,]
  
  if (grepl(".csv$",countfile)) {
    data<-read.csv(countfile,header=T,row.names=1,as.is=T,check.names=FALSE)
  } else {
    data<-read.delim(countfile,header=T,row.names=1,as.is=T,check.names=FALSE)
  }
  
  if(transformTable){
    data<-t(data)
  }
  
  data<-data[,colnames(data) != "Feature_length"]
  colClass<-sapply(data, class)
  countNotNumIndex<-which(colClass!="numeric" & colClass!="integer")
  if (length(countNotNumIndex)==0) {
    index<-1;
    indecies<-c()
  } else {
    index<-max(countNotNumIndex)+1
    indecies<-c(1:(index-1))
  }
  
  countData<-data[,c(index:ncol(data))]
  countData[is.na(countData)] <- 0
  countData<-round(countData)
  
  if(addCountOne){
    countData<-countData+1
  }
  
  comparisonNames=comparisons$ComparisonName
  
  pairedspearman<-list()
  
  newVarInData<-setdiff(colnames(data),colnames(dataAllOut))
  if (length(newVarInData)>0) {
    dataAllOut<-align(dataAllOut,data[,newVarInData,drop=FALSE])
  }
  resultAllOutVar<-c("baseMean","log2FoldChange","pvalue","padj")
  
  comparison_index = 1
  
  for(comparison_index in c(1:nrow(comparisons))){
    comparisonName=comparisons$ComparisonName[comparison_index]
    comparisonTitle=comparisons$ComparisonTitle[comparison_index]
    titles<-c(titles, comparisonTitle)
    cat(comparisonName, " ", comparisonTitle, "\n")
    designFile=comparisons$ConditionFile[comparison_index]
    gnames=unlist(comparisons[comparison_index, c("ReferenceGroupName", "SampleGroupName")])
    
    designData<-read.table(designFile, sep="\t", header=T)
    designData$Condition<-factor(designData$Condition, levels=gnames)
    
    if(ncol(designData) >= 3){
      cat("Data with covariances!\n")
    }else{
      cat("Data without covariances!\n")
    }
    if (any(colnames(designData)=="Paired")) {
      ispaired<-TRUE
      cat("Paired Data!\n")
    }else{
      ispaired<-FALSE
      cat("Not Paired Data!\n")
    }
    temp<-apply(designData,2,function(x) length(unique(x)))
    if (any(temp==1)) {
      cat(paste0("Factors with only 1 level in design matrix: ",colnames(designData)[which(temp==1)],"\n"))
      cat("They will be removed")
      cat("\n")
      designData<-designData[,which(temp!=1)]
    }
    temp<-apply(designData[,-1,drop=F],2,rank)
    if (length(unique(rowSums(temp)))==1 | identical(temp[,1],temp[,-1])) {
      cat(paste0("The model matrix is not full rank, so the model cannot be fit as specified"))
      cat("\n")
      cat("Only Condition variable will be kept.")
      cat("\n")
      designData<-designData[,which(colnames(designData)%in% c("Sample","Condition"))]
    }
    
    comparisonData<-countData[,colnames(countData) %in% as.character(designData$Sample),drop=F]
    if(ncol(comparisonData) != nrow(designData)){
      message=paste0("Data not matched, there are ", nrow(designData), " samples in design file ", designFile, " but ", ncol(comparisonData), " samples in data ")
      warning(message)
      writeLines(message,paste0(comparisonName,".error"))
      next
    }
    comparisonData<-comparisonData[,as.character(designData$Sample)]
    
    prefix<-paste0(comparisonName, suffix)
    if(top25only){
      ranks=apply(comparisonData, 2, function(x){
        y=x[x > 0]
        q=quantile(y)
        return(x>=q[4])
      })
      
      select=apply(ranks, 1, function(x){
        any(x)
      })
      
      comparisonData=comparisonData[select,]
    }
    
    if(detectedInBothGroup){
      conds<-unique(designData$Condition)
      data1<-comparisonData[, colnames(comparisonData) %in% designData$Sample[designData$Condition==conds[1]],drop=FALSE]
      data2<-comparisonData[, colnames(comparisonData) %in% designData$Sample[designData$Condition==conds[2]],drop=FALSE]
      med1<-apply(data1, 1, median) > zeroCount
      med2<-apply(data2, 1, median) > zeroCount
      med<-med1 & med2
      
      comparisonData<-comparisonData[med,]
    }
    
    if(performWilcox){
      #quantile and wilcox
      quantileData=normalize.quantiles(data.matrix(comparisonData))
      colnames(quantileData)=colnames(comparisonData)
      rownames(quantileData)=rownames(comparisonData)
      write.csv(quantileData, file=paste0(prefix, "_quantile.csv"), row.names = T)
      
      data1<-quantileData[, colnames(quantileData) %in% designData$Sample[designData$Condition==conds[1]],drop=FALSE]
      data2<-quantileData[, colnames(quantileData) %in% designData$Sample[designData$Condition==conds[2]],drop=FALSE]
      
      diffData=data.frame(quantileData)
      diffData$pvalues=unlist(lapply(c(1:nrow(data1)), function(index){
        d1=data1[index,]
        d2=data2[index,]
        test=wilcox.test(d1,d2)
        test$p.value
      }))
      diffData$log2MedianFoldChange=unlist(lapply(c(1:nrow(data1)), function(index){
        d1=data1[index,]
        d2=data2[index,]
        log2(median(d2) / median(d1))
      }))
      diffData$log2MeanFoldChange=unlist(lapply(c(1:nrow(data1)), function(index){
        d1=data1[index,]
        d2=data2[index,]
        log2(mean(d2) / mean(d1))
      }))
      diffData=diffData[order(diffData$pvalues),]
      write.csv(diffData, file=paste0(prefix, "_quantile_wilcox.csv"), row.names = T)
      
      filterData=diffData[diffData$pvalues<=pvalue & abs(diffData$log2MedianFoldChange) > log2(foldChange),]
      write.csv(filterData, file=paste0(prefix, "_quantile_wilcox_sig.csv"), row.names = T)
    }
    
    if(minMedianInGroup > 0){
      conds<-unique(designData$Condition)
      data1<-comparisonData[, colnames(comparisonData) %in% designData$Sample[designData$Condition==conds[1]],drop=FALSE]
      data2<-comparisonData[, colnames(comparisonData) %in% designData$Sample[designData$Condition==conds[2]],drop=FALSE]
      med1<-apply(data1, 1, median) >= minMedianInGroup
      med2<-apply(data2, 1, median) >= minMedianInGroup
      med<-med1 | med2
      
      comparisonData<-comparisonData[med,]
      
      cat(nrow(comparisonData), " genes with minimum median count in group larger or equals than ", minMedianInGroup, "\n")
    }
    
    if (nrow(comparisonData)==0) {
      message=paste0("Error: 0 Genes can be used in DESeq2 analysis in comparison ",comparisonName," \n")
      warning(message)
      writeLines(message,paste0(comparisonName,".error"))
      next;
    }
    
    validComparisons<-c(validComparisons, comparisonName)
    
    if(ispaired){
      pairedSamples = unique(designData$Paired)
      
      spcorr<-unlist(lapply(c(1:length(pairedSamples)), function(x){
        samples<-designData$Sample[designData$Paired==pairedSamples[x]]
        cor(comparisonData[,samples[1]],comparisonData[,samples[2]],method="spearman")
      }))
      
      
      sptable<-data.frame(Name=pairedSamples, Spcorr=spcorr)
      write.csv(sptable, file=paste0(prefix, "_Spearman.csv"), row.names=FALSE)

      dir.create("details", showWarnings = FALSE)
      
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
    
    if(ispaired){
      colnames(comparisonData)<-unlist(lapply(c(1:ncol(comparisonData)), function(i){paste0(designData$Paired[i], "_", colnames(comparisonData)[i])}))
    }
    rownames(designData)<-colnames(comparisonData)
    conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[designData$Condition]))
    
    write.csv(comparisonData, file=paste0(prefix, ".csv"))
    
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
    
    width=max(4000, ncol(rldmatrix) * 40 + 1000)
    height=max(3000, ncol(rldmatrix) * 40)
    png(filename=paste0(prefix, "_DESeq2-log2-density-individual.png"), width=width, height=height, res=300)
    g<-ggplot(rsdata) + geom_density(aes(x=log2Count, colour=Sample)) + facet_wrap(~Sample, scales = "free") + xlab("DESeq2 log2 transformed count")
    print(g)
    dev.off()
    
    
    #varianceStabilizingTransformation
    
    dds<-myEstimateSizeFactors(dds)
    
    fitType<-"parametric"
    if(nrow(comparisonData) < 5){
      fitType<-"mean"
    }
    while(1){
      #varianceStabilizingTransformation
      vsdres<-try(vsd <- varianceStabilizingTransformation(dds, blind=TRUE,fitType=fitType))
      if(class(vsdres) == "try-error"){
        if(grepl("every gene contains at least one zero", vsdres[1])){
          removed<-removed+1
          keptNumber<-length(zeronumbers) - percent10 * removed
          keptSample<-zeronumbers[1:keptNumber]
          excludedSample<-zeronumbers[(keptNumber+1):length(zeronumbers)]
          
          comparisonData<-comparisonData[, colnames(comparisonData) %in% keptSample]
          designData<-designData[rownames(designData) %in% keptSample,]
          dds=DESeqDataSetFromMatrix(countData = comparisonData,
                                     colData = designData,
                                     design = ~1)
          
          colnames(dds)<-colnames(comparisonData)
        } else if (grepl("newsplit: out of vertex space", vsdres[1]) | fitType != "mean") {
          message=paste0("Warning: varianceStabilizingTransformation function can't run. fitType was set to mean to try again")
          warning(message)
          fitType<-"mean"
          writeLines(message,paste0(comparisonName,".error"))
        } else {
          message=paste0(paste0("Error: varianceStabilizingTransformation function can't run. ", vsdres))
          writeLines(message,paste0(comparisonName,".error"))
          stop(message)
        }
      }else if(all(is.na(assay(vsd)))){
        fitType<-"mean"
      } else{
        conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[designData$Condition]))
        break
      }
    }
    
    if(nrow(comparisonData) > 1){
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
    }
    
    #different expression analysis
    designFormula=as.formula(paste0("~",paste0(c(colnames(designData)[-c(1:2)],"Condition"),collapse="+")))
    
    cat(paste0("", designFormula), "\n")
    
    dds=DESeqDataSetFromMatrix(countData = comparisonData,
                               colData = designData,
                               design = designFormula)
    
    dds<-myEstimateSizeFactors(dds)
    
    ddsres<-try(dds <- DESeq(dds,fitType=fitType))
    if(class(ddsres) == "try-error"){
      if( grepl("One can instead use the gene-wise estimates as final estimates", ddsres[1])){
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds<-nbinomWaldTest(dds)
      }else if(grepl("newsplit: out of vertex space", ddsres[1])){
        dds <- DESeq(dds,fitType="mean")
      }else{
        stop(paste0("DESeq2 failed: ", ddsres[1]))
      }
    }
    res<-results(dds,cooksCutoff=FALSE)
    
    cat("DESeq2 finished.\n")
    
    if (useRawPvalue==1) {
      select<-(!is.na(res$pvalue)) & (res$pvalue<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))
    } else {
      select<-(!is.na(res$padj)) & (res$padj<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))
    }
    
    if(length(indecies) > 0){
      inddata<-data[rownames(comparisonData),indecies,drop=F]
      tbb<-cbind(inddata, as.data.frame(comparisonData), res)
    }else{
      tbb<-cbind(as.data.frame(comparisonData), res)
    }

    tbb$FoldChange<-2^tbb$log2FoldChange
    tbbselect<-tbb[select,,drop=F]
    tbbAllOut<-as.data.frame(tbb[,resultAllOutVar,drop=F])
    tbbAllOut$Significant<-select
    colnames(tbbAllOut)<-paste0(colnames(tbbAllOut)," (",comparisonName,")")
    resultAllOut<-cbind(as.data.frame(resultAllOut)[row.names(dataAllOut),],as.matrix(tbbAllOut[row.names(dataAllOut),]))
    row.names(resultAllOut)<-row.names(dataAllOut)
    
    tbb<-tbb[order(tbb$pvalue),,drop=F]
    write.csv(as.data.frame(tbb),paste0(prefix, "_DESeq2.csv"))
    
    tbbselect<-tbbselect[order(tbbselect$pvalue),,drop=F]
    sigFile=paste0(prefix, "_DESeq2_sig.csv")
    sigTable<-as.data.frame(tbbselect)
    write.csv(sigTable,sigFile)
    
    allSigNameList[[comparisonName]]<-row.names(sigTable)
    allSigDirectionList[[comparisonName]]<-sign(sigTable$log2FoldChange)
    if(nrow(sigTable) > 0){
      sigTable$comparisonName<-comparisonName
      
      if (("Feature_gene_name" %in% colnames(sigTable)) & (!("Feature_gene_name" %in% sigTableAllVar))){
        sigTableAllVar<-c("Feature_gene_name", sigTableAllVar)
      }
      
      sigTableAll<-rbind(sigTableAll,sigTable[,c("comparisonName",sigTableAllVar),drop=FALSE],make.row.names=FALSE)
      sigTableAllGene<-c(sigTableAllGene,row.names(sigTable))
    }
    
    if("Feature_gene_name" %in% colnames(tbb)){
      write.table(tbb[,c("Feature_gene_name", "stat"),drop=F],paste0(prefix, "_DESeq2_GSEA.rnk"),row.names=F,col.names=F,sep="\t", quote=F)
      if(exportSignificantGeneName){
        write.table(tbbselect[,c("Feature_gene_name"),drop=F], paste0(prefix, "_DESeq2_sig_genename.txt"),row.names=F,col.names=F,sep="\t", quote=F)
      }
    }else{
      write.table(tbb[,c("stat"),drop=F],paste0(prefix, "_DESeq2_GSEA.rnk"),row.names=T,col.names=F,sep="\t", quote=F)
      if(exportSignificantGeneName){
        write.table(data.frame(name=rownames(tbbselect)), paste0(prefix, "_DESeq2_sig_genename.txt"),row.names=F,col.names=F,sep="\t", quote=F)
      }
    }    
    
    if(showDEGeneCluster){
      siggenes<-rownames(rldmatrix) %in% rownames(tbbselect)
      
      nonDEmatrix<-rldmatrix[!siggenes,,drop=F]
      DEmatrix<-rldmatrix[siggenes,,drop=F]
      
      drawPCA(paste0(prefix,"_geneDE"),DEmatrix , showLabelInPCA, designData, conditionColors)
      drawHCA(paste0(prefix,"_geneDE"),DEmatrix , ispaired, designData, conditionColors, gnames)
      
      drawPCA(paste0(prefix,"_geneNotDE"), nonDEmatrix, showLabelInPCA, designData, conditionColors)
      drawHCA(paste0(prefix,"_geneNotDE"), nonDEmatrix, ispaired, designData, conditionColors, gnames)
      
      #drawHCA(paste0(prefix,"_gene500NotDE"), nonDEmatrix[1:min(500, nrow(nonDEmatrix)),,drop=F], ispaired, designData, conditionColors, gnames)
    }
    
    #Top 25 Significant genes barplot
    sigDiffNumber<-nrow(tbbselect)
    if (sigDiffNumber>0) {
      if (sigDiffNumber>25) {
        print(paste0("More than 25 genes were significant. Only the top 25 genes will be used in barplot"))
        diffResultSig<-tbbselect[order(tbbselect$pvalue)[1:25],]
      } else {
        diffResultSig<-tbbselect
      }
      if("Feature_gene_name" %in% colnames(diffResultSig)){
        diffResultSig$Name<-as.character(diffResultSig$Feature_gene_name)
      }else{
        diffResultSig$Name<-sapply(strsplit(row.names(diffResultSig),";"),function(x) x[1])
      }
      diffResultSig$Name <- factor(diffResultSig$Name, levels=diffResultSig$Name[order(diffResultSig$log2FoldChange)])
      diffResultSig<-as.data.frame(diffResultSig)
      
      png(filename=paste0(prefix, "_DESeq2_sig_barplot.png"), width=3000, height=3000, res=300)
      #   pdf(paste0(prefix,"_DESeq2_sig_barplot.pdf"))
      p<-ggplot(diffResultSig,aes(x=Name,y=log2FoldChange,order=log2FoldChange))+geom_bar(stat="identity")+
        coord_flip()+
        #     geom_abline(slope=0,intercept=1,colour="red",linetype = 2)+
        scale_y_continuous(name=bquote(log[2]~Fold~Change))+
        theme(axis.text = element_text(colour = "black"))
      print(p)
      dev.off()
    } else {
      print(paste0("No gene with adjusted p value less than ",pvalue," and fold change larger than ",foldChange))
    }
    
    #volcano plot
    changeColours<-c(grey="grey",blue="blue",red="red")
    diffResult<-as.data.frame(tbb)
    diffResult$log10BaseMean<-log10(diffResult$baseMean)
    diffResult$colour<-"grey"
    if (useRawPvalue==1) {
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    } else {
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    }
    
    png(filename=paste0(prefix, "_DESeq2_volcanoPlot.png"), width=3000, height=3000, res=300)
    #  pdf(paste0(prefix,"_DESeq2_volcanoPlot.pdf"))
    if (useRawPvalue==1) {
      p<-ggplot(diffResult,aes(x=log2FoldChange,y=pvalue))+
        scale_y_continuous(trans=reverselog_trans(10),name=bquote(p~value))
    } else {
      p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
        scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))
    }
    p<-p+geom_point(aes(size=log10BaseMean,colour=colour))+
      scale_color_manual(values=changeColours,guide = FALSE)+
      scale_x_continuous(name=bquote(log[2]~Fold~Change))+
      geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
      geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
      guides(size=guide_legend(title=bquote(log[10]~Base~Mean)))+
      theme_bw()+
      scale_size(range = c(3, 7))+
      theme(axis.text = element_text(colour = "black",size=30),
            axis.title = element_text(size=30),
            legend.text= element_text(size=30),
            legend.title= element_text(size=30))
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
}

allprefix=basename(inputfile)
#Venn for all significant genes
#Output all significant genes table
if(!is.null(sigTableAll)){
  sigTableAll<-cbind(Gene=sigTableAllGene,sigTableAll)
  write.csv(sigTableAll,paste0(allprefix, "_DESeq2_allSig.csv"),row.names=FALSE)
  
  #Do venn if length between 2-5
  if (length(allSigNameList)>=2 & length(allSigNameList)<=5) {
    venn.diagram1<-function (x, filename, height = 3000, width = 3000, resolution = 500, 
                             units = "px", compression = "lzw", na = "stop", main = NULL, 
                             sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
                             main.fontfamily = "serif", main.col = "black", main.cex = 1, 
                             main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                             sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
                             sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE,
                             fill=NA,
                             ...) 
    {
      if (is.na(fill[1])) {
        if (length(x)==5) {
          fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")
        } else if (length(x)==4) {
          fill = c("dodgerblue", "goldenrod1",  "seagreen3", "orchid3")
        } else if (length(x)==3) {
          fill = c("dodgerblue", "goldenrod1", "seagreen3")
        } else if (length(x)==2) {
          fill = c("dodgerblue", "goldenrod1")
        }
      }
      if (force.unique) {
        for (i in 1:length(x)) {
          x[[i]] <- unique(x[[i]])
        }
      }
      if ("none" == na) {
        x <- x
      }
      else if ("stop" == na) {
        for (i in 1:length(x)) {
          if (any(is.na(x[[i]]))) {
            stop("NAs in dataset", call. = FALSE)
          }
        }
      }
      else if ("remove" == na) {
        for (i in 1:length(x)) {
          x[[i]] <- x[[i]][!is.na(x[[i]])]
        }
      }
      else {
        stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
      }
      if (0 == length(x) | length(x) > 5) {
        stop("Incorrect number of elements.", call. = FALSE)
      }
      if (1 == length(x)) {
        list.names <- category.names
        if (is.null(list.names)) {
          list.names <- ""
        }
        grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
                                                   category = list.names, ind = FALSE,fill=fill, ...)
      }
      else if (2 == length(x)) {
        grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
                                                     area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
                                                                                                           x[[2]])), category = category.names, ind = FALSE, 
                                                     fill=fill,
                                                     ...)
      }
      else if (3 == length(x)) {
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        list.names <- category.names
        nab <- intersect(A, B)
        nbc <- intersect(B, C)
        nac <- intersect(A, C)
        nabc <- intersect(nab, C)
        grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
                                                   area2 = length(B), area3 = length(C), n12 = length(nab), 
                                                   n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
                                                   category = list.names, ind = FALSE, list.order = 1:3, 
                                                   fill=fill,
                                                   ...)
      }
      else if (4 == length(x)) {
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        D <- x[[4]]
        list.names <- category.names
        n12 <- intersect(A, B)
        n13 <- intersect(A, C)
        n14 <- intersect(A, D)
        n23 <- intersect(B, C)
        n24 <- intersect(B, D)
        n34 <- intersect(C, D)
        n123 <- intersect(n12, C)
        n124 <- intersect(n12, D)
        n134 <- intersect(n13, D)
        n234 <- intersect(n23, D)
        n1234 <- intersect(n123, D)
        grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
                                                 area2 = length(B), area3 = length(C), area4 = length(D), 
                                                 n12 = length(n12), n13 = length(n13), n14 = length(n14), 
                                                 n23 = length(n23), n24 = length(n24), n34 = length(n34), 
                                                 n123 = length(n123), n124 = length(n124), n134 = length(n134), 
                                                 n234 = length(n234), n1234 = length(n1234), category = list.names, 
                                                 ind = FALSE, fill=fill,...)
      }
      else if (5 == length(x)) {
        A <- x[[1]]
        B <- x[[2]]
        C <- x[[3]]
        D <- x[[4]]
        E <- x[[5]]
        list.names <- category.names
        n12 <- intersect(A, B)
        n13 <- intersect(A, C)
        n14 <- intersect(A, D)
        n15 <- intersect(A, E)
        n23 <- intersect(B, C)
        n24 <- intersect(B, D)
        n25 <- intersect(B, E)
        n34 <- intersect(C, D)
        n35 <- intersect(C, E)
        n45 <- intersect(D, E)
        n123 <- intersect(n12, C)
        n124 <- intersect(n12, D)
        n125 <- intersect(n12, E)
        n134 <- intersect(n13, D)
        n135 <- intersect(n13, E)
        n145 <- intersect(n14, E)
        n234 <- intersect(n23, D)
        n235 <- intersect(n23, E)
        n245 <- intersect(n24, E)
        n345 <- intersect(n34, E)
        n1234 <- intersect(n123, D)
        n1235 <- intersect(n123, E)
        n1245 <- intersect(n124, E)
        n1345 <- intersect(n134, E)
        n2345 <- intersect(n234, E)
        n12345 <- intersect(n1234, E)
        grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), 
                                                      area2 = length(B), area3 = length(C), area4 = length(D), 
                                                      area5 = length(E), n12 = length(n12), n13 = length(n13), 
                                                      n14 = length(n14), n15 = length(n15), n23 = length(n23), 
                                                      n24 = length(n24), n25 = length(n25), n34 = length(n34), 
                                                      n35 = length(n35), n45 = length(n45), n123 = length(n123), 
                                                      n124 = length(n124), n125 = length(n125), n134 = length(n134), 
                                                      n135 = length(n135), n145 = length(n145), n234 = length(n234), 
                                                      n235 = length(n235), n245 = length(n245), n345 = length(n345), 
                                                      n1234 = length(n1234), n1235 = length(n1235), n1245 = length(n1245), 
                                                      n1345 = length(n1345), n2345 = length(n2345), n12345 = length(n12345), 
                                                      category = list.names, ind = FALSE,fill=fill, ...)
      }
      else {
        stop("Invalid size of input object")
      }
      if (!is.null(sub)) {
        grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, 
                               fontface = sub.fontface, fontfamily = sub.fontfamily, 
                               col = sub.col, cex = sub.cex)
      }
      if (!is.null(main)) {
        grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, 
                               fontface = main.fontface, fontfamily = main.fontfamily, 
                               col = main.col, cex = main.cex)
      }
      grid.newpage()
      grid.draw(grob.list)
      return(1)
      # return(grob.list)
    }
    makeColors<-function(n,colorNames="Set1") {
      maxN<-brewer.pal.info[colorNames,"maxcolors"]
      if (n<=maxN) {
        colors<-brewer.pal(n, colorNames)
        if (length(colors)>n) {
          colors<-colors[1:n]
        }
      } else {
        colors<-colorRampPalette(brewer.pal(maxN, colorNames))(n)
      }
      return(colors)
    }
    colors<-makeColors(length(allSigNameList))
    png(paste0(allprefix,"_significantVenn.png"),res=300,height=2000,width=2000)
    venn.diagram1(allSigNameList,cex=2,cat.cex=2,cat.col=colors,fill=colors)
    dev.off()
  }
  #Do heatmap significant genes if length larger or equal than 2
  if (length(allSigNameList)>=2) {
    temp<-cbind(unlist(allSigNameList),unlist(allSigDirectionList))
    colnames(temp)<-c("Gene","Direction")
    temp<-cbind(temp,comparisonName=rep(names(allSigNameList),sapply(allSigNameList,length)))
    temp<-data.frame(temp)
    dataForFigure<-temp
    #geting dataForFigure order in figure
    temp$Direction<-as.integer(as.character(temp$Direction))
    temp<-acast(temp, Gene~comparisonName ,value.var="Direction")
    temp<-temp[do.call(order, data.frame(temp)),]
    maxNameChr<-max(nchar(row.names(temp)))
    if (maxNameChr>70) {
      row.names(temp)<-substr(row.names(temp),0,70)
      dataForFigure$Gene<-substr(dataForFigure$Gene,0,70)
      warning(paste0("The gene names were too long (",maxNameChr,"). Only first 70 letters were kept."))
    }
    dataForFigure$Gene<-factor(dataForFigure$Gene,levels=row.names(temp))
    
    width=max(2500, 60 * length(unique(dataForFigure$comparisonName)))
    height=max(2000, 40 * length(unique(dataForFigure$Gene)))
    png(paste0(allprefix, "_significantHeatmap.png"),res=300,height=height,width=width)
    g<-ggplot(dataForFigure, aes(comparisonName, Gene))+
      geom_tile(aes(fill=Direction), color="white") +
      scale_fill_manual(values=c("light green", "red")) +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0.5, face="bold"),
            axis.text.y = element_text(size=textSize, face="bold")) +
      coord_equal()
    print(g)
    dev.off()
  }
}

if (! is.null(resultAllOut)) {
  #write a file with all information
  resultAllOut<-cbind(dataAllOut,resultAllOut[row.names(dataAllOut),])
  write.csv(resultAllOut,paste0(allprefix, "_DESeq2.csv"))
  
  if(length(validComparisons) > 1 ){
    #volcano plot for all comparisons
    temp<-resultAllOut[,-(1:ncol(dataAllOut))]
    diffResult<-NULL
    diffResultVar<-unique(sapply(strsplit(colnames(temp)," "),function(x) x[1]))
    for (i in 1:(length(validComparisons))) {
      temp1<-temp[,(i*length(diffResultVar)-(length(diffResultVar)-1)):(i*length(diffResultVar))]
      colnames(temp1)<-diffResultVar
      temp1$Comparison<-validComparisons[i]
      if (is.null(diffResult)) {
        diffResult<-temp1
      } else {
        diffResult<-rbind(diffResult,temp1)
      }
    }
    changeColours<-c(grey="grey",blue="blue",red="red")
    diffResult$log10BaseMean<-log10(diffResult$baseMean)
    diffResult$Comparison<-allTitles[diffResult$Comparison]
    diffResult$Comparison<-factor(diffResult$Comparison,levels=unique(diffResult$Comparison))
    diffResult$colour<-"grey"
    if (useRawPvalue==1) {
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    } else {
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    }
    
    width<-max(2000,2000*length(allComparisons))
    png(filename=paste0(allprefix, "_DESeq2_volcanoPlot.png"), width=width, height=2000, res=300)
    #  pdf(paste0(prefix,"_DESeq2_volcanoPlot.pdf"))
    if (useRawPvalue==1) {
      p<-ggplot(diffResult,aes(x=log2FoldChange,y=pvalue))+
        scale_y_continuous(trans=reverselog_trans(10),name=bquote(p~value))
    } else {
      p<-ggplot(diffResult,aes(x=log2FoldChange,y=padj))+
        scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))
    }
    p<-p+geom_point(aes(size=log10BaseMean,colour=colour))+
      scale_color_manual(values=changeColours,guide = FALSE)+
      scale_x_continuous(name=bquote(log[2]~Fold~Change))+
      geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
      geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
      guides(size=guide_legend(title=bquote(log[10]~Base~Mean)))+
      theme_bw()+
      scale_size(range = c(3, 7))+
      facet_grid(. ~ Comparison)+
      theme(axis.text = element_text(colour = "black",size=30),
            axis.title = element_text(size=30),
            legend.text= element_text(size=30),
            legend.title= element_text(size=30),
            strip.text.x = element_text(size = 30))
    print(p)
    dev.off()
	
	#output a summary table with numbers of gisnificant changed genes
	sigGeneSummaryTable<-t(table(diffResult[,"Significant"],diffResult[,"Comparison"]))
	sigGeneSummaryTable<-data.frame(Comparison=row.names(sigGeneSummaryTable),GeneInComparison=rowSums(sigGeneSummaryTable),NotSignificant=sigGeneSummaryTable[,1],Significant=sigGeneSummaryTable[,2])
	write.csv(sigGeneSummaryTable,paste0(allprefix, "_DESeq2_sigGeneSummary.csv"),row.names=FALSE)
	
  }  
}

