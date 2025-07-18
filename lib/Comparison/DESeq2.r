rootdir<-"/nobackup/shah_lab/shengq2/20250113_APD_Wora_Project/20250113_12615_hg38/deseq2_genetable/result"
inputfile<-"APD_Wora.define" 

pvalue<-0.1
useRawPvalue<-0
foldChange<-1.2
minMedianInGroup<-5
  
detectedInBothGroup<-0
showLabelInPCA<-0
showDEGeneCluster<-1
addCountOne<-0
usePearsonInHCA<-1
top25only<-0
performWilcox<-0
textSize<-10
transformTable<-0
exportSignificantGeneName<-1
thread<-8

independentFiltering<-TRUE

outputPdf<-FALSE;outputPng<-TRUE;outputTIFF<-FALSE;showVolcanoLegend<-FALSE;top25cvInHCA<-FALSE;
cooksCutoff<-TRUE
#predefined_condition_end

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

if(useRawPvalue){
  alpha<-0.1
  suffix=paste0(suffix, "_pvalue", pvalue)
}else{
  alpha<-pvalue
  suffix=paste0(suffix, "_fdr", pvalue)
}

zeroCount=0
if(addCountOne){
  zeroCount=1
  minMedianInGroup=minMedianInGroup+1
}

if(!exists("idIndex")){
  idIndex<-1
}

if(!exists("outputPdf")){
  outputPdf<-FALSE
}

if(!exists("outputPng") | !outputPdf ){
  outputPng<-TRUE
}

if(!exists("outputTIFF")){
  outputTIFF<-FALSE
}

if(!exists("filterBaseMean")){
  filterBaseMean<-0
}

if(!exists("filterBaseMeanValue")){
	filterBaseMeanValue<-30
}

outputFormat<-c()
if(outputPdf){
  outputFormat<-c("PDF")
}
if(outputPng){
  outputFormat<-c(outputFormat, "PNG")
}
if(outputTIFF){
  outputFormat<-c(outputFormat, "TIFF")
}

if(!exists("countSep")){
  countSep="\t"
}

if(!exists("usePearsonInHCA")){
  usePearsonInHCA=0
}

if(!exists("exportSignificantGeneName")){
  exportSignificantGeneName<-1
}

if(exists("libraryFile")){
  if (libraryKey != 'None'){
    if (grepl(".csv$", libraryFile)){
      librarySize_tbl<-read.csv(libraryFile, row.names=1,check.names=FALSE)
      if(libraryKey == "sizeFactor"){
        librarySize<-unlist(librarySize_tbl[,libraryKey,drop=T])
        names(librarySize)<-rownames(librarySize_tbl)
      }else{
        librarySize<-unlist(librarySize_tbl[libraryKey,,drop=T])
      }
      cat("Using ", libraryKey, " in " , libraryFile , " as library size. \n")
    }else{
      librarySize<-read.table(libraryFile, row.names=1,check.names=FALSE,header=T,stringsAsFactor=F)
      librarySize<-unlist(librarySize[,libraryKey,drop=T])
      cat("Using ", libraryKey, " in " , libraryFile , " as library size. \n")
    }
  }
}

if(!exists("thread")){
  thread<-1
}

if(!exists("showVolcanoLegend")){
  showVolcanoLegend<-0
}

if(!exists("cooksCutoff")){
  cooksCutoff<-"DEFAULT"
}

cat("cooksCutoff=",cooksCutoff,"\n")

library("DESeq2")
library("lattice")
library("ggplot2")
library("grid")
library("scales")
library("reshape2")
library("VennDiagram")
library("RColorBrewer")
library("BiocParallel")
library("ggrepel")
library("stringr")
library("data.table")
library("cowplot")
suppressPackageStartupMessages(library("ComplexHeatmap"))

setwd(rootdir)  

source('countTableVisFunctions.R')

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

theme_bw3 <- function (axis.x.rotate=F) { 
  is_ggplot2_newver = packageVersion("ggplot2") >= "3.4.0"

  if(is_ggplot2_newver){
    eline = element_line(colour = "black", linewidth = 0.5)
  }else{
    eline = element_line(colour = "black", size = 0.5)
  }

  result = theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      plot.title = element_text(hjust = 0.5),
      axis.line = eline
    )

  if (axis.x.rotate){
    result = result + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  return(result)
}

# Convert a byte-compiled function to an interpreted-code function 
unByteCode <- function(fun) {
  FUN <- eval(parse(text=deparse(fun)))
  environment(FUN) <- environment(fun)
  FUN
}

# Replace function definition inside of a locked environment **HACK** 
assignEdgewise <- function(name, env, value) {
  unlockBinding(name, env=env)
  assign( name, envir=env, value=value)
  lockBinding(name, env=env)
  invisible(value)
}

# Replace byte-compiled function in a locked environment with an interpreted-code
# function
unByteCodeAssign <- function(fun) {
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

openPlot<-function(filePrefix, format, width_inch, height_inch, figureName){
  fileName<-paste0(filePrefix, ".", tolower(format))
  if(format == "PDF"){
    pdf(fileName, width=width_inch, height=height_inch, useDingbats=FALSE)
  }else if(format == "TIFF"){
    tiff(filename=fileName, width=width_inch, height=height_inch, units="in", res=300)
  }else {
    png(filename=fileName, width=width_inch, height=height_inch, units="in", res=300)
  }
  cat("saving", figureName, "to ", fileName, "\n")
}

drawPlot<-function(filePrefix, outputFormat, width_inch, height_inch, p, figureName){
  for(format in outputFormat){  
    openPlot(filePrefix, format, width_inch, height_inch, figureName)  
    print(p)
    dev.off()
  }
}

drawHCA<-function(prefix, rldselect, ispaired, designData, groupColors, gnames, outputFormat, legend_label_gp, column_names_gp){
  mat_scaled = t(scale(t(rldselect)))

  ha=HeatmapAnnotation( Group=designData$Condition,
                        col=list(Group=groupColors),
                        annotation_legend_param = list(Group = list(ncol = 1, 
                                                                    title = "Group", 
                                                                    title_position = "topleft"),
                                                                    title_gp=legend_label_gp, 
                                                                    labels_gp=legend_label_gp))

  if(usePearsonInHCA){
    clustering_distance_columns="pearson"
  }else{
    clustering_distance_columns="euclidean"
  }

  draw_heatmap_png( filepath=paste0(prefix, "_DESeq2-vsd-heatmap.png"), 
                    htdata=mat_scaled, 
                    name="Z-Score", 
                    show_row_names=FALSE, 
                    show_column_names=TRUE,
                    save_rds=FALSE,
                    clustering_distance_columns=clustering_distance_columns,
                    top_annotation=ha,
                    show_row_dend=FALSE,
                    column_names_gp = column_names_gp,
                    legend_gp = legend_label_gp)
}

myEstimateSizeFactors<-function(dds){
  if(exists("librarySize")){
    cat("Estimate size factor based on library size\n")
    curLibrarySize<-librarySize[colnames(dds)]
    if(libraryKey == "sizeFactor"){
      curSizeFactor<- curLibrarySize
    }else{
      #based on DESeq2 introduction
      curSizeFactor<- curLibrarySize / exp(mean(log(curLibrarySize)))
    }
    sizeFactors(dds)<-curSizeFactor
  }else{
    cat("Estimate size factor based on reads\n")
    sfres<-try(dds<-estimateSizeFactors(dds))
    if (class(sfres) == "try-error") {
      library(edgeR)
      countNum<-counts(dds)
      y<-calcNormFactors(countNum, methold="TMM")
      cs<-colSums(countNum)
      cs<-cs / median(cs)
      sf<-y * cs
      sizeFactors(dds)<-sf
    }
  }
  return(dds)
}

###########################
#end function
###########################
# 
# ###################################################################
# #change comparisons_data, need to be removed before adding to pipeline
# comparisons_data=rbind(comparisons_data,comparisons_data)
# comparisons_data[3:4,1]=c("Control_placenta_vs_Heart","Diabetic_placenta_vs_Heart")
# comparisons_data[3:4,6]=c("Control_placenta_vs_Heart","Diabetic_placenta_vs_Heart")
# comparisons_data[,1]=paste0("Test_",comparisons_data[,1])
# comparisons_data[,3]="/scratch/cqs/zhaos/RolandaLister/20200907_RolandaLister4363_4369_RnaSeq/pipeline/deseq2_proteincoding_genetable/result/test.design"
# comparisons_data[,6]=paste0("Test_",comparisons_data[,6])
# comparisons_data$designFormula="~Tissue + Condition+Tissue:Condition"
# comparisons_data$contrast=c("Condition_Diabetic_vs_Control",paste0("Condition_Diabetic_vs_Control",";","Tissueplacenta.ConditionDiabetic"),
#                             "Tissue_placenta_vs_Heart",paste0("Tissue_placenta_vs_Heart",";","Tissueplacenta.ConditionDiabetic"))
# #comparisons_data=comparisons_data[1,,drop=FALSE]
# #end change comparisons_data
# ###################################################################

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

feature_name_regex=NA
feature_name_filter=NA
n_first=-1
enhanced_volcano_red_blue_only=FALSE
title_in_volcano=TRUE
caption_in_volcano=TRUE

heatmap_add_width_inch=2
heatmap_add_height_inch=0

heatmap_legend_label_fontsize=18
heatmap_column_name_fontsize=18

de_biotype=NA
selectLab=NULL

DE_combatseq = 0

if(file.exists("fileList1.txt")){
  options_table = read.table("fileList1.txt", sep="\t")
  myoptions = split(options_table$V1, options_table$V2)
  feature_name_regex = myoptions$feature_name_regex
  feature_name_filter = myoptions$feature_name_filter
  enhanced_volcano_red_blue_only = myoptions$enhanced_volcano_red_blue_only == "1"
  title_in_volcano = myoptions$title_in_volcano == "1"
  caption_in_volcano = myoptions$caption_in_volcano == "1"
  if("n_first" %in% names(myoptions)){
    n_first = as.numeric(myoptions$n_first)
  }

  heatmap_add_width_inch=to_numeric(myoptions$heatmap_add_width_inch, heatmap_add_width_inch)
  heatmap_add_height_inch=to_numeric(myoptions$heatmap_add_height_inch, heatmap_add_height_inch)

  heatmap_legend_label_fontsize=to_numeric(myoptions$heatmap_legend_label_fontsize, heatmap_legend_label_fontsize)
  heatmap_column_name_fontsize=to_numeric(myoptions$heatmap_column_name_fontsize, heatmap_column_name_fontsize)

  de_biotype = myoptions$de_biotype
  if(!is.na(de_biotype)){
    if(de_biotype == ""){
      de_biotype = NA
    }
  }

  if(!is_one(myoptions$showLabelInVolcano)){
    selectLab=""
  }

  DE_combatseq = is_one(myoptions$DE_combatseq)
}

if(!is.na(de_biotype)){
  de_biotype_name=gsub(" ", "_", de_biotype)
  cat("DE biotype: ", de_biotype, "\n")
}
legend_label_gp = gpar(fontsize = heatmap_legend_label_fontsize, fontface = "bold")
column_names_gp = gpar(fontsize = heatmap_column_name_fontsize, fontface = "bold")

countfile_index = 1
titles<-NULL
validComparisons<-c()
for(countfile_index in c(1:length(countfiles))){
  countfile = countfiles[countfile_index]
  comparisons = comparisons_data[comparisons_data$CountFile == countfile,]
  
  if(n_first != -1){
    data<-data.frame(fread(countfile, nrows=n_first), row.names=idIndex,check.names=FALSE)
  }else{
    if (grepl(".csv$",countfile)) {
      data<-read.csv(countfile,header=T,row.names=idIndex,as.is=T,check.names=FALSE)
    } else {
      data<-read.delim(countfile,header=T,row.names=idIndex,as.is=T,check.names=FALSE, sep=countSep)
    }
  }
  
  if(transformTable){
    data<-t(data)
  }

  if(!is.na(feature_name_filter)){
    if(!is.null(feature_name_filter)){
      if(feature_name_filter != ""){
        if("Feature_gene_name" %in% colnames(data)){
          cat("Filter feature name by ", feature_name_filter, "\n")
          discarded=data[!grepl(feature_name_filter, data$Feature_gene_name, perl=TRUE),]
          cat("  Discarded ", nrow(discarded), " features\n")
          data=data[grepl(feature_name_filter, data$Feature_gene_name, perl=TRUE),]
        }
      }
    }
  }

  if(!is.na(feature_name_regex)){
    if(!is.null(feature_name_regex)){
      if(feature_name_regex != ""){
        rownames(data) = str_match(rownames(data), feature_name_regex)[,2]
      }
    }
  }

  if(is.na(de_biotype)){
    de_features=c()
  }else{
    de_features=rownames(data)[data$Feature_gene_biotype == de_biotype]
  }

  data<-data[,colnames(data) != "Feature_length"]
  colClass<-sapply(data, class)
  countNotNumIndex<-which((colClass!="numeric" & colClass!="integer") | grepl("Gene_Id", colnames(data)))
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

    if ("pairOnlyCovariant" %in% colnames(comparisons)) {
      pairOnlyCovariant=comparisons$pairOnlyCovariant[comparison_index]
      if (is.na(pairOnlyCovariant) || (pairOnlyCovariant=="")) {
        pairOnlyCovariant=NULL
      }
    }else{
      pairOnlyCovariant=NULL
    }

    if ("designFormula" %in% colnames(comparisons)) {
      designFormula=comparisons$designFormula[comparison_index]
      print(paste0("designFormula = ", designFormula, "\n"))
      if (is.na(designFormula) || (designFormula=="")) {
        designFormula=NULL
      } else {
        designFormula=as.formula(designFormula)
      }
    } else {
      designFormula=NULL
    }

    if ("contrast" %in% colnames(comparisons)) {
      contrast=comparisons$contrast[comparison_index]
      if (is.na(contrast) || (contrast=="")) {
        contrast=NULL
      } else {
        contrast=list(strsplit(contrast,";")[[1]])
      }
    } else {
      contrast=NULL
    }

    if ("collapse_by" %in% colnames(comparisons)) {
      collapse_by=comparisons$collapse_by[comparison_index]
      if (is.na(collapse_by) || (collapse_by=="")) {
        collapse_by=NULL
      }
    } else {
      collapse_by=NULL
    }
    
    titles<-c(titles, comparisonTitle)
    cat(comparisonName, " ", comparisonTitle, "\n")
    designFile=comparisons$ConditionFile[comparison_index]
    #comment here as has many group names
    gnames=unlist(comparisons[comparison_index, c("ReferenceGroupName", "SampleGroupName")])
    #gnames=as.character(unique(designData$Condition))
    groupColors=c("blue", "red")
    names(groupColors)<-gnames
    
    designData<-read.table(designFile, sep="\t", header=T)
    designData$Condition<-factor(designData$Condition, levels=gnames)

    if(!is.null(pairOnlyCovariant)){
      if(!any(colnames(designData) == pairOnlyCovariant)){
        stop(paste0("Cannot find pairOnlyCovariant ", pairOnlyCovariant, " in ",designFile))
      }
      tbl = table(designData[,pairOnlyCovariant])
      tbl = tbl[tbl == 2]
      designData=designData[designData[,pairOnlyCovariant] %in% names(tbl),,drop=F]
    }
    
    missedSamples<-as.character(designData$Sample)[!(as.character(designData$Sample) %in% colnames(countData))]
    if(length(missedSamples) > 0){
      message=paste0("There are missed sample defined in design file but not in real data: ", missedSamples)
      warning(message)
      writeLines(message,paste0(comparisonName,".error"))
      next
    }

    comparisonData<-countData[,colnames(countData) %in% as.character(designData$Sample),drop=F]
    if(ncol(comparisonData) != nrow(designData)){
      message=paste0("Data not matched, there are ", nrow(designData), " samples in design file ", designFile, " but ", ncol(comparisonData), " samples in data ")
      warning(message)
      writeLines(message,paste0(comparisonName,".error"))
      next
    }
    comparisonData<-comparisonData[,as.character(designData$Sample)]

    if(!is.null(collapse_by)){
      dds=DESeqDataSetFromMatrix(countData = comparisonData,
                                colData = designData,
                                design = ~1) 
      dds=collapseReplicates(dds, designData[,collapse_by], designData$Sample)   
      designData<-designData[!duplicated(designData[,collapse_by]),]
      designData$Sample<-designData[,collapse_by]
      designData<-designData[,colnames(designData) != collapse_by]
      comparisonData<-counts(dds)[,designData$Sample]
      rm(dds)
    }

    if("batch" %in% colnames(designData) & DE_combatseq){
      library(sva)
      cat("CombatSeq: ", comparisonName, "\n")
      if(ncol(designData) > 3){
        cat("Data with covariances!\n")
        covar_mod=designData |>
          dplyr::select(-Sample, -Condition, -batch)
      }else{
        cat("Data without covariances!\n")
        covar_mod=NULL
      }
      comparisonData <- ComBat_seq(comparisonData, batch=designData$batch, group=designData$Condition, covar_mod=covar_mod)
      write.csv(comparisonData, file=paste0(comparisonName, "_combatseq.csv"))
      cat("Remove batch from the design matrix\n")
      designData<-designData |> dplyr::select(-batch)
    }
    
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
      
      geneNumBeforeFilter=nrow(comparisonData)
      comparisonData<-comparisonData[med,]
      
      cat(nrow(comparisonData), " genes with minimum median count in group larger or equals than ", minMedianInGroup, ". ",geneNumBeforeFilter-nrow(comparisonData)," genes removed\n")
    }
    
    if (nrow(comparisonData)<=1) {
      message=paste0("Error: Only ", nrow(comparisonData), " Genes can be used in DESeq2 analysis in comparison ",comparisonName,", ignored. \n")
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
    
    write.csv(comparisonData, file=paste0(prefix, ".csv"))
    
    #some basic graph
    dds=DESeqDataSetFromMatrix(countData = comparisonData,
                               colData = designData,
                               design = ~1)
    
    colnames(dds)<-colnames(comparisonData)
		dds<-myEstimateSizeFactors(dds)

    if(filterBaseMean){
      cat(paste0("filter by basemean: ", filterBaseMeanValue, "\n"))
      baseMeans = rowMeans(counts(dds, normalized=TRUE))
      write.csv(baseMeans, file=paste0(prefix, ".basemean.csv"))
      dds<-dds[baseMeans > filterBaseMeanValue,]
      comparisonData=comparisonData[baseMeans > filterBaseMeanValue,]
    }

    rld_normed<-log2(counts(dds,normalized=TRUE) + 1)
    write.csv(format(rld_normed, digits=3), paste0(prefix, "_DESeq2-log2-normalized-counts.csv"))

    #draw density graph
    log2counts<-as.matrix(log2(counts(dds,normalized=FALSE) + 1))

    if(is.na(de_biotype)){
      draw_density_plot(log2counts, paste0(prefix, "_DESeq2-log2-density"), outputFormat)
    }else{
      biotype_log2counts=log2counts[rownames(log2counts) %in% de_features,]
      draw_density_plot(biotype_log2counts, paste0(prefix, "_", de_biotype_name, "_DESeq2-log2-density"), outputFormat)
    }

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
        break
      }
    }
    
    assayvsd<-assay(vsd)
    write.csv(format(assayvsd, digits=3), file=paste0(prefix, "_DESeq2-vsd.csv"))
    
    rldmatrix=as.matrix(assayvsd)

    if(nrow(comparisonData) > 1){
      #draw pca graph
      if(is.na(de_biotype)){
        drawPCA(file_prefix=paste0(prefix,"_geneAll_DESeq2-vsd-pca"),
                rldmatrix=rldmatrix, 
                showLabelInPCA=showLabelInPCA, 
                groups=designData$Condition, 
                groupColors=groupColors, 
                outputFormat=outputFormat, 
                scalePCs=TRUE)

        if(exists("top25cvInHCA") && top25cvInHCA){
          rv<-rowVars(rldmatrix)
          countHT<-rldmatrix[rv>=quantile(rv)[4],]
          drawHCA(prefix=paste0(prefix,"_geneTop25variance"), 
                  rldselect=countHT, 
                  ispaired=ispaired, 
                  designData=designData, 
                  groupColors=groupColors, 
                  gnames=ganems, 
                  outputFormat=outputFormat,
                  legend_label_gp=legend_label_gp,
                  column_names_gp=column_names_gp)
        }else{
          drawHCA(prefix=paste0(prefix,"_geneAll"), 
                  rldselect=rldmatrix, 
                  ispaired=ispaired, 
                  designData=designData, 
                  groupColors=groupColors, 
                  gnames=ganems, 
                  outputFormat=outputFormat,
                  legend_label_gp=legend_label_gp,
                  column_names_gp=column_names_gp)
        }
      }else{
        biotype_rldmatrix=rldmatrix[rownames(rldmatrix) %in% de_features,]
        drawPCA(file_prefix=paste0(prefix,"_", de_biotype_name, "_DESeq2-vsd-pca"),
                rldmatrix=rldmatrix, 
                showLabelInPCA=showLabelInPCA, 
                groups=designData$Condition, 
                groupColors=groupColors, 
                outputFormat=outputFormat, 
                scalePCs=TRUE)

        drawHCA(prefix=paste0(prefix,"_", de_biotype_name), 
                rldselect=biotype_rldmatrix, 
                ispaired=ispaired, 
                designData=designData, 
                groupColors=groupColors, 
                gnames=ganems, 
                outputFormat=outputFormat,
                legend_label_gp=legend_label_gp,
                column_names_gp=column_names_gp)
      }
    }

    #different expression analysis
    if (is.null(designFormula)) {
      designFormula=as.formula(paste0("~",paste0(c(colnames(designData)[-c(1:2)],"Condition"),collapse="+")))
    }
    
    cat(paste0("", designFormula), "\n")
    
    dds=DESeqDataSetFromMatrix(countData = comparisonData,
                               colData = designData,
                               design = designFormula)
    
    dds<-myEstimateSizeFactors(dds)

    if(!is.na(de_biotype)){
      #we need to use all genes for size factor estimation but only care about certain biotype
      dds = dds[rownames(dds) %in% de_features,]
      de_prefix = paste0(prefix, "_", de_biotype_name)
      cat("Only use", nrow(dds), de_biotype, "genes for DE analysis\n")
    }else{
      de_prefix = prefix
    }

    bpparam<-MulticoreParam(thread)
#    parallel<-ifelse(thread <= 1, FALSE, TRUE)
    parallel=FALSE
    
    ddsres<-try(dds <- DESeq(dds,fitType=fitType, parallel=parallel, BPPARAM=bpparam))
    
    if(class(ddsres) == "try-error"){
      if( grepl("One can instead use the gene-wise estimates as final estimates", ddsres[1])){
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds<-nbinomWaldTest(dds)
      }else if(grepl("newsplit: out of vertex space", ddsres[1])){
        dds <- DESeq(dds,fitType="mean", parallel=parallel, BPPARAM=bpparam)
      }else{
        stop(paste0("DESeq2 failed: ", ddsres[1]))
      }
    }
    
    cat("Using independentFiltering =", independentFiltering, "\n")
    if(cooksCutoff == "TRUE" || cooksCutoff == "DEFAULT"){
      cat("Using default cooksCutoff\n")
      if (!is.null(contrast)) {
        res<-results(dds, alpha=alpha, parallel=parallel, BPPARAM=bpparam, independentFiltering=independentFiltering, contrast=contrast) 
      } else {
        res<-results(dds, alpha=alpha, parallel=parallel, BPPARAM=bpparam, independentFiltering=independentFiltering)
      }
    }else{
      cat("Using cooksCutoff =", cooksCutoff, "\n")
      if (!is.null(contrast)) {
        res<-results(dds, cooksCutoff=cooksCutoff, alpha=alpha, parallel=parallel, BPPARAM=bpparam, independentFiltering=independentFiltering, contrast=contrast) 
      } else {
        res<-results(dds, cooksCutoff=cooksCutoff, alpha=alpha, parallel=parallel, BPPARAM=bpparam, independentFiltering=independentFiltering)
      }
    }

    res$FoldChange<-2^res$log2FoldChange

    baseMeanPerLvl <- sapply( levels(dds$Condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$Condition == lvl,drop=FALSE] ) )
    colnames(baseMeanPerLvl)<-paste0("baseMean_", colnames(baseMeanPerLvl))
    res<-cbind(res, baseMeanPerLvl)

    cat("DESeq2 finished.\n")
    
    if (useRawPvalue==1) {
      select<-(!is.na(res$pvalue)) & (res$pvalue<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))
    } else {
      select<-(!is.na(res$padj)) & (res$padj<pvalue) & ((res$log2FoldChange >= log2(foldChange)) | (res$log2FoldChange <= -log2(foldChange)))
    }
    
    if(length(indecies) > 0){
      inddata<-data[rownames(res),indecies,drop=F]
      tbb<-cbind(inddata, as.data.frame(comparisonData[rownames(res),,drop=FALSE]), res)
    }else{
      tbb<-cbind(as.data.frame(comparisonData[rownames(res),,drop=FALSE]), res)
    }
    
    tbbselect<-tbb[select,,drop=F]
    tbbselect$DE_status=ifelse(tbbselect$log2FoldChange > 0, "UP", "DOWN")
    tbb$DE_status="NS"
    tbb[rownames(tbbselect),"DE_status"]<-tbbselect$DE_status

    tbbAllOut<-as.data.frame(tbb[,resultAllOutVar,drop=F])
    tbbAllOut$Significant<-select
    colnames(tbbAllOut)<-paste0(colnames(tbbAllOut)," (",comparisonName,")")

    if(0){
      dl = list(resultAllOut=resultAllOut, dataAllOut=dataAllOut, tbbAllOut=tbbAllOut)
      saveRDS(dl, paste0(prefix, ".debug.rds"))

      d1=readRDS('AtRisk1_2_vs_AtRisk0_min5_fdr0.05.debug.rds')
      dl=readRDS("AtRisk1_vs_AtRisk0_min5_fdr0.05.debug.rds")
      resultAllOut<-dl$resultAllOut
      dataAllOut<-dl$dataAllOut
      tbbAllOut<-dl$tbbAllOut
    }
    
    if(is.null(resultAllOut)){
      resultAllOut<-merge(dataAllOut |> tibble::rownames_to_column("Feature"),
                          tbbAllOut |> tibble::rownames_to_column("Feature"),
                          by="Feature", 
                          all=TRUE)
    }else{
      resultAllOut<-merge(resultAllOut, 
                          tbbAllOut |> tibble::rownames_to_column("Feature"), 
                          by="Feature", 
                          all=TRUE)
    }
    
    tbb<-tbb[order(tbb$pvalue),,drop=F]
    write.csv(as.data.frame(tbb),paste0(de_prefix, "_DESeq2.csv"))
    
    tbbselect<-tbbselect[order(tbbselect$pvalue),,drop=F]
    sigFile=paste0(de_prefix, "_DESeq2_sig.csv")
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
    
    geneNameField = NULL
    
    lowColNames = tolower(colnames(tbb))
    for(name in c("Feature_gene_name", "Gene.Symbol", "Gene_Symbol", "Gene Symbol", "gene_name")){
      lowName = tolower(name)
      if(lowName %in% lowColNames){
        geneNameField=colnames(tbb)[match(lowName, lowColNames)]
        break
      }
    }
    
    if(!is.null(geneNameField)){
      write.table(tbb[,c(geneNameField, "stat"),drop=F],paste0(de_prefix, "_DESeq2_GSEA.rnk"),row.names=F,col.names=F,sep="\t", quote=F)
      if(exportSignificantGeneName){
        write.table(tbbselect[,c(geneNameField),drop=F], paste0(de_prefix, "_DESeq2_sig_genename.txt"),row.names=F,col.names=F,sep="\t", quote=F)
      }
    }else{
      write.table(tbb[,c("stat"),drop=F],paste0(de_prefix, "_DESeq2_GSEA.rnk"),row.names=T,col.names=F,sep="\t", quote=F)
      if(exportSignificantGeneName){
        write.table(data.frame(name=rownames(tbbselect)), paste0(de_prefix, "_DESeq2_sig_genename.txt"),row.names=F,col.names=F,sep="\t", quote=F)
      }
    }    
    
    if(showDEGeneCluster){
      siggenes<-rownames(rldmatrix) %in% rownames(tbbselect)
      
      nonDEmatrix<-rldmatrix[!siggenes,,drop=F]
      DEmatrix<-rldmatrix[siggenes,,drop=F]
      
      drawPCA(file_prefix=paste0(de_prefix,"_geneDE_DESeq2-vsd-pca"),
              rldmatrix=DEmatrix, 
              showLabelInPCA=showLabelInPCA, 
              groups=designData$Condition, 
              groupColors=groupColors, 
              outputFormat=outputFormat, 
              scalePCs=TRUE)

      drawHCA(paste0(de_prefix,"_geneDE"),
              DEmatrix, 
              ispaired, 
              designData, 
              groupColors=groupColors, 
              gnames, 
              outputFormat,
              legend_label_gp=legend_label_gp,
              column_names_gp=column_names_gp)
      
      drawPCA(file_prefix=paste0(de_prefix,"_geneNotDE_DESeq2-vsd-pca"),
              rldmatrix=nonDEmatrix, 
              showLabelInPCA=showLabelInPCA, 
              groups=designData$Condition, 
              groupColors=groupColors, 
              outputFormat=outputFormat, 
              scalePCs=TRUE)

      drawHCA(paste0(de_prefix,"_geneNotDE"), 
              nonDEmatrix, 
              ispaired, 
              designData, 
              groupColors=groupColors, 
              gnames, 
              outputFormat,
              legend_label_gp=legend_label_gp,
              column_names_gp=column_names_gp)
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
      if(!is.null(geneNameField)){
        diffResultSig$Name<-as.character(diffResultSig[,geneNameField])
      }else{
        diffResultSig$Name<-sapply(strsplit(row.names(diffResultSig),";"),function(x) x[1])
      }
      if (any(duplicated(diffResultSig$Name))) {
        whichIndex<-which(duplicated(diffResultSig$Name))
        diffResultSig$Name[whichIndex]<-paste0(row.names(diffResultSig)[whichIndex], ":", diffResultSig$Name[whichIndex])
      }
      diffResultSig$Name <- factor(diffResultSig$Name, levels=diffResultSig$Name[order(diffResultSig$log2FoldChange)])
      diffResultSig<-as.data.frame(diffResultSig)
      
      p<-ggplot(diffResultSig,aes(x=Name,y=log2FoldChange,order=log2FoldChange))+
        geom_bar(stat="identity")+
        coord_flip()+
        #     geom_abline(slope=0,intercept=1,colour="red",linetype = 2)+
        scale_y_continuous(name=bquote(log[2](fold~change)))+
        theme_bw3() +
        theme(axis.text = element_text(colour = "black"),
              axis.title.x = element_blank())
      
      filePrefix<-paste0(de_prefix,"_DESeq2_sig_barplot")
      drawPlot(filePrefix, outputFormat, 6, 6, p, "PCA")
    } else {
      print(paste0("No gene with adjusted p value less than ",pvalue," and fold change larger than ",foldChange))
    }
    
    #volcano plot
    changeColours<-c(grey="grey",blue="blue",red="red")
    diffResult<-as.data.frame(tbb)
    diffResult$log10BaseMean<-log10(diffResult$baseMean)
    diffResult$colour<-"grey"
    if (useRawPvalue==1) {
      diffResult<-subset(diffResult, !is.na(pvalue))
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    } else {
      diffResult<-subset(diffResult, !is.na(padj))
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    }

    diffResult$colour=factor(diffResult$colour,levels=c("grey","blue","red"))
    diffResult$log10pvalue=-log10(diffResult$pvalue)
    diffResult<-diffResult[order(diffResult$colour),]
    diffResult$colour<-as.character(diffResult$colour)

    write.csv(diffResult, file=paste0(de_prefix, "_DESeq2_volcanoPlot.csv"))
    
    xname=bquote(log[2](fold~change))
    yname=bquote(-log[10](p~value))

    p<-ggplot(diffResult,aes(x=log2FoldChange,y=log10pvalue))+
      scale_x_continuous(name=xname,breaks=pretty_breaks(n=4))+
      scale_y_continuous(name=yname,breaks=pretty_breaks(n=4))+
      geom_point(aes(size=log10BaseMean,colour=colour))+
      scale_color_manual(values=changeColours,guide = FALSE)+
      geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
      geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
      guides(size=guide_legend(title=bquote(log[10](BaseMean))))+
      theme_bw()+
      scale_size(range = c(3, 7))+
      theme(axis.text = element_text(colour = "black",size=30),
            axis.title = element_text(size=30),
            legend.text= element_text(size=30),
            legend.title= element_text(size=30))
    
    if(!showVolcanoLegend){
      p<-p+ theme(legend.position = "none")
      width_inch=7
    }else{
      width_inch=10
    }
    
    filePrefix<-paste0(de_prefix,"_DESeq2_volcanoPlot")
    drawPlot(filePrefix, outputFormat, width_inch, 7, p, "Volcano")

    if(require("EnhancedVolcano", quietly=TRUE)){
      if(!("Feature_gene_name" %in% colnames(diffResult))){
        diffResult$Feature_gene_name=rownames(diffResult)
      }

      if(enhanced_volcano_red_blue_only){
        keyvals.colour=diffResult$colour
        names(keyvals.colour)[keyvals.colour == 'red'] <- 'Up'
        names(keyvals.colour)[keyvals.colour == 'grey'] <- 'NS'
        names(keyvals.colour)[keyvals.colour == 'blue'] <- 'Down'
      }else{
        keyvals.colour=NULL
      }

      if(title_in_volcano){
        title=gsub("_vs_", " vs ", comparisonTitle)
      }else{
        title=NULL
      }

      if(caption_in_volcano){
        caption=paste0("total = ", nrow(diffResult), " variables")
      }else{
        caption=NULL
      }

      diffResult$Short_name<-gsub(";.+","", diffResult$Feature_gene_name)
      if(packageVersion("EnhancedVolcano") == '1.8.0'){
        if(useRawPvalue == 1){
          yvar="pvalue"
        }else{
          yname=bquote(-log[10](adjusted~p~value))
          yvar="padj"
        }

        p<-EnhancedVolcano(diffResult,
            lab = diffResult$Short_name,
            selectLab = selectLab,
            x = 'log2FoldChange',
            y = yvar,
            pCutoff = pvalue,
            FCcutoff = log2(foldChange),
            pointSize = 3.0,
            labSize = 6.0,
            colCustom = keyvals.colour,
            colAlpha = 1,
            title = title,
            subtitle = NULL,
            caption = caption)
      }else{
        if(useRawPvalue == 1){
          pCutoffCol="pvalue"
        }else{
          pCutoffCol="padj"
        }
        p<-EnhancedVolcano(diffResult,
            lab = diffResult$Short_name,
            selectLab = selectLab,
            x = 'log2FoldChange',
            y = 'pvalue',
            pCutoff = pvalue,
            pCutoffCol = pCutoffCol,
            FCcutoff = log2(foldChange),
            pointSize = 3.0,
            labSize = 6.0,
            colCustom = keyvals.colour,
            colAlpha = 1,
            title = title,
            subtitle = NULL,
            caption = caption) 
      }
      p=p + 
        ylab(yname) + 
        xlab(xname) +
        theme(plot.title = element_text(hjust = 0.5))
      filePrefix<-paste0(de_prefix,"_DESeq2_volcanoEnhanced")
      drawPlot(filePrefix, outputFormat, 7, 7, p, "Volcano")
    }
  }
  
  if(length(pairedspearman) > 0){
    filePrefix<-paste0(de_prefix, "_", ifelse(minMedianInGroup > 0, paste0("spearman_min", minMedianInGroup), "spearman"))
    fwidth<-max(6, 1 * length(pairedspearman))
    for(format in outputFormat){
      openPlot(filePrefix, format, fwidth, 6, "Spearman correlation")
      boxplot(pairedspearman)
      dev.off()
    }
  }
}

allprefix=paste0(basename(inputfile), suffix)
if(!is.na(de_biotype)){
  allprefix=paste0(allprefix, "_", de_biotype_name)
  resultAllOut=resultAllOut |> dplyr::filter(Feature_gene_biotype == de_biotype)
}

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
    for(format in outputFormat){ 
      filePrefix<-paste0(allprefix,"_significantVenn")
      openPlot(filePrefix, format, 7, 7, "Venn")
      venn.diagram1(allSigNameList,cex=2,cat.cex=2,cat.col=colors,fill=colors)
      dev.off()
    }
  }
  # #Do heatmap significant genes if length larger or equal than 2
  # if (length(allSigNameList)>=2) {
  #   temp<-cbind(unlist(allSigNameList),unlist(allSigDirectionList))
  #   colnames(temp)<-c("Gene","Direction")
  #   temp<-cbind(temp,comparisonName=rep(names(allSigNameList),sapply(allSigNameList,length)))
  #   temp<-data.frame(temp)
  #   dataForFigure<-temp
  #   #geting dataForFigure order in figure
  #   temp$Direction<-as.integer(as.character(temp$Direction))
  #   temp<-acast(temp, Gene~comparisonName ,value.var="Direction")
  #   temp<-temp[do.call(order, data.frame(temp)),]
  #   maxNameChr<-max(nchar(row.names(temp)))
  #   if (maxNameChr>70) {
  #     tmpNames<-substr(row.names(temp),0,70)
  #     if(length(tmpNames) == length(unique(tmpNames))){
  #       row.names(temp)<-tmpNames
  #       dataForFigure$Gene<-substr(dataForFigure$Gene,0,70)
  #       warning(paste0("The gene names were too long (",maxNameChr,"). Only first 70 letters were kept."))
  #     }
  #   }
  #   dataForFigure$Gene<-factor(dataForFigure$Gene,levels=row.names(temp))
    
  #   g<-ggplot(dataForFigure, aes(comparisonName, Gene))+
  #     geom_tile(aes(fill=Direction), color="white") +
  #     scale_fill_manual(values=c("light green", "red")) +
  #     theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0.5, face="bold"),
  #           axis.text.y = element_text(size=textSize, face="bold")) +
  #     coord_equal()
    
  #   width=min(max(8, 3 * length(unique(dataForFigure$comparisonName))),20)
  #   height=min(max(7, 40 * length(unique(dataForFigure$Gene))),100)
  #   filePrefix<-paste0(allprefix,"_significantHeatmap")
  #   drawPlot(filePrefix, outputFormat, 7, 7, g, "Significant Heatmap")
  # }
}

if (! is.null(resultAllOut)) {
  #write a file with all information
  resultAllOut<-resultAllOut |> tibble::column_to_rownames("Feature")
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
      diffResult<-subset(diffResult, !is.na(pvalue))
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$pvalue<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    } else {
      diffResult<-subset(diffResult, !is.na(padj))
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange>=log2(foldChange))]<-"red"
      diffResult$colour[which(diffResult$padj<=pvalue & diffResult$log2FoldChange<=-log2(foldChange))]<-"blue"
    }

    diffResult$colour=factor(diffResult$colour,levels=c("grey","blue","red"))
    diffResult$log10pvalue=-log10(diffResult$pvalue)

    filePrefix<-paste0(allprefix,"_DESeq2_volcanoPlot")
    saveRDS(diffResult, paste0(filePrefix, ".rds"))
    
    pair_font_size=25
    p<-ggplot(diffResult,aes(x=log2FoldChange,y=log10pvalue))+
      geom_point(aes(colour=colour))+
      scale_color_manual(values=changeColours,guide = FALSE)+
      scale_x_continuous(name=bquote(log[2](fold~change)))+
      scale_y_continuous(name=bquote(-log[10](p~value))) +
      geom_hline(yintercept = 1,colour="grey",linetype = "dotted")+
      geom_vline(xintercept = 0,colour="grey",linetype = "dotted")+
      theme_bw()+
      facet_grid(. ~ Comparison)+
      theme(axis.text = element_text(colour = "black",size=25),
            axis.title = element_text(size=25),
            legend.text= element_text(size=25),
            legend.title= element_text(size=25),
            strip.text.x = element_text(size=pair_font_size),
            strip.background=element_rect(fill="white"))

    max_pair=max(get_text_width(unique(diffResult$Comparison), font_family = "", font_size = pair_font_size)) + 0.1

    pwidth<-max(12, length(unique(diffResult$Comparison))*max_pair+1)
    drawPlot(filePrefix, outputFormat, pwidth, 7, p, "Volcano")
    
    #output a summary table with numbers of gisnificant changed genes
    sigGeneSummaryTable<-t(table(diffResult[,"Significant"],diffResult[,"Comparison"]))
    
    notSigIndex<-match("FALSE", colnames(sigGeneSummaryTable))
    if(is.na(notSigIndex)){
      notSignificant=0
    }else{
      notSignificant=sigGeneSummaryTable[,notSigIndex]
    }
    
    sigIndex<-match("TRUE", colnames(sigGeneSummaryTable))
    if(is.na(sigIndex)){
      significant=0
    }else{
      significant=sigGeneSummaryTable[,sigIndex]
    }
    
    dSigGeneSummaryTable<-data.frame(Comparison=row.names(sigGeneSummaryTable),GeneInComparison=rowSums(sigGeneSummaryTable),NotSignificant=notSignificant,Significant=significant)
    write.csv(dSigGeneSummaryTable,paste0(allprefix, "_DESeq2_sigGeneSummary.csv"),row.names=FALSE)
  }  
}

#export session information
writeLines(capture.output(sessionInfo()), paste0(basename(inputfile),".DESeq2.SessionInfo.txt"))
deseq2version<-paste0("DESeq2,v", packageVersion("DESeq2"))

