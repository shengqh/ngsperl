
library(limma)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(scales)

########################################################
# Functions
########################################################

makeDesignTable <- function(rawDataTable, factorTable,factorInDesign=NULL,addIntercept=1) {
  if (is.null(factorInDesign)) {
    factorInDesign=as.character(unique(factorTable[, 2]))
  }
  designTable <- matrix(0, nrow = ncol(rawDataTable), ncol = length(factorInDesign))
  row.names(designTable) <- colnames(rawDataTable)
  colnames(designTable) <- factorInDesign
  
  designTable[as.matrix(factorTable)] <- 1
  if (addIntercept) {
    designTable <- cbind(Intercept = 1, designTable)
  }
  designTable=designTable[which(rowSums(designTable)>0),]
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


#maybe need put these into a seprate file as these are analysis functions not "report" functions
library("VennDiagram")
venn.diagram1<-function (x, count=NULL,filename, height = 3000, width = 3000, resolution = 500, 
                         units = "px", compression = "lzw", na = "stop", main = NULL, 
                         sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
                         main.fontfamily = "serif", main.col = "black", main.cex = 1, 
                         main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                         sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
                         sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE,
                         fill=NA,
                         ...) 
{
  if (is.null(count)) {
    countFun<-function(x) length(x)
  } else {
    countFun<-function(x) sum(count[x])
  }
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
    grob.list <- VennDiagram::draw.single.venn(area = countFun(x[[1]]), 
                                               category = list.names, ind = FALSE,fill=fill, ...)
  }
  else if (2 == length(x)) {
    grob.list <- VennDiagram::draw.pairwise.venn(area1 = countFun(x[[1]]), 
                                                 area2 = countFun(x[[2]]), cross.area = countFun(intersect(x[[1]], 
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
    grob.list <- VennDiagram::draw.triple.venn(area1 = countFun(A), 
                                               area2 = countFun(B), area3 = countFun(C), n12 = countFun(nab), 
                                               n23 = countFun(nbc), n13 = countFun(nac), n123 = countFun(nabc), 
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
    grob.list <- VennDiagram::draw.quad.venn(area1 = countFun(A), 
                                             area2 = countFun(B), area3 = countFun(C), area4 = countFun(D), 
                                             n12 = countFun(n12), n13 = countFun(n13), n14 = countFun(n14), 
                                             n23 = countFun(n23), n24 = countFun(n24), n34 = countFun(n34), 
                                             n123 = countFun(n123), n124 = countFun(n124), n134 = countFun(n134), 
                                             n234 = countFun(n234), n1234 = countFun(n1234), category = list.names, 
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
    grob.list <- VennDiagram::draw.quintuple.venn(area1 = countFun(A), 
                                                  area2 = countFun(B), area3 = countFun(C), area4 = countFun(D), 
                                                  area5 = countFun(E), n12 = countFun(n12), n13 = countFun(n13), 
                                                  n14 = countFun(n14), n15 = countFun(n15), n23 = countFun(n23), 
                                                  n24 = countFun(n24), n25 = countFun(n25), n34 = countFun(n34), 
                                                  n35 = countFun(n35), n45 = countFun(n45), n123 = countFun(n123), 
                                                  n124 = countFun(n124), n125 = countFun(n125), n134 = countFun(n134), 
                                                  n135 = countFun(n135), n145 = countFun(n145), n234 = countFun(n234), 
                                                  n235 = countFun(n235), n245 = countFun(n245), n345 = countFun(n345), 
                                                  n1234 = countFun(n1234), n1235 = countFun(n1235), n1245 = countFun(n1245), 
                                                  n1345 = countFun(n1345), n2345 = countFun(n2345), n12345 = countFun(n12345), 
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
  #	return(grob.list)
}

saveFigure=function(file,type="png",width=2000,height=2000,res=300) {
  png(file,width=width,height=height,res=res)
}


########################################################
# Parameters
########################################################

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
rawData <- read.csv(rawDataFile, header = T, row.names = 1, as.is = T, check.names=F)

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
saveFigure(paste0(outFile,".RawBoxplot.png"))
boxplot(rawDataTable, las = 2)
dev.off()

dataForPlot=reshape2::melt(rawDataTable,variable.name="Sample")
p=ggplot(dataForPlot,aes(x=value,colour=Sample))+geom_density()
saveFigure(paste0(outFile,".RawDensity.png"))
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
  saveFigure(paste0(outFile,".NormlizationBoxplot.png"))
  boxplot(rawDataTable, las = 2)
  dev.off()
  
  dataForPlot=reshape2::melt(rawDataTable,variable.name="Sample")
  p=ggplot(dataForPlot,aes(x=value,colour=Sample))+geom_density()
  saveFigure(paste0(outFile,".NormlizationDensity.png"))
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

saveFigure(paste0(outFile,".Heatmap.png"))
print(ht)
dev.off()

##################################
# save overall design table
##################################
if (!exists("noFactorSampleAsIntercept") || noFactorSampleAsIntercept==1) { #use samples without factor define as Intercept in design
  designTableOverall=makeDesignTable(rawDataTable, factorTable,addIntercept=1)
} else {
  designTableOverall=makeDesignTable(rawDataTable, factorTable,addIntercept=0)
}
write.csv(designTableOverall,paste0(outFile, ".OverallDesignTable.csv"))


##################################
# save differential criteria
##################################
#diffCriteriaTable=cbind(Criteria=c("Perform Normlization","NA sample propotion in Group","Use Raw p value","p value","Fold change"),Value=c(performNormlization,notNaPropotion,useRawPvalue,pvalue,foldChange))
diffCriteriaTable=data.frame(performNormlization,notNaPropotion,useRawPvalue,pvalue,foldChange)
write.csv(diffCriteriaTable,paste0(outFile, ".DiffCriteria.csv"),row.names=FALSE)


##################################
# For each comprasion, perform limma
##################################
naCount=apply(rawDataTable,1,function(x) length(which(is.na(x))))
rawDataTableForDiffDetection=rawDataTable[which(naCount/ncol(rawDataTable)<=1-notNaPropotion),]

DiffResultsOutTableAll=list() #record all results when needed
DiffResultsOutDiffGenesAll=list() #record all diff genes when needed
comprasionsTable=NULL #record comprasions information
for (ComparisonOne in unique(factorInComparisons[,3])) {
  print(paste0("Working in ",ComparisonOne))
  #  prefix=paste0(outFile,"_",ComparisonOne)
  prefix=paste0(ComparisonOne)
  
  factorInDesign=unique(factorInComparisons[which(factorInComparisons[,3]==ComparisonOne & factorInComparisons[,2]!="contrast"),1])
  if (!exists("noFactorSampleAsIntercept") || noFactorSampleAsIntercept==1) { #use samples without factor define as Intercept in design
    rawDataTableForDiffDetectionOne=rawDataTableForDiffDetection
    designTable <- makeDesignTable(rawDataTableForDiffDetectionOne, factorTable,factorInDesign=factorInDesign)
  } else { #remove samples without factor define
    rawDataTableForDiffDetectionOne=rawDataTableForDiffDetection[,unique(factorTable[which(factorTable[,2] %in% factorInDesign),1])]
    designTable <- makeDesignTable(rawDataTableForDiffDetectionOne, factorTable,factorInDesign=factorInDesign,addIntercept=0)
  }
  
  
  if (any(factorInComparisons[,3]==ComparisonOne & factorInComparisons[,2]=="contrast")) {
    contrastsText=factorInComparisons[which(factorInComparisons[,3]==ComparisonOne & factorInComparisons[,2]=="contrast")[1],1]
  } else {
    if (exists("noFactorSampleAsIntercept") && noFactorSampleAsIntercept==0) { #remove samples without factor define and no Intercept in design, contrast maybe needed
      warning(paste0("No contrast defined and no Intercept in design, contrast maybe needed to get correct comparison"))
    }
    contrastsText=NULL
  }
  #record comprasions information in comprasionsTable
  comprasionsTable=rbind(comprasionsTable,c(ComparisonOne,paste(factorInDesign,collapse=";"),ifelse(is.null(contrastsText),factorInDesign[1],contrastsText)))
  
  fit3 <- performLimma(rawDataTableForDiffDetectionOne, designTable,contrastsText=contrastsText)
  DiffResults <- topTable(fit3, coef = 1,number=99999)
  
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
  tbb=tbb %>% mutate(Significant=case_when(tbb[,log2FcCol]>0 ~ "SignificantUp",tbb[,log2FcCol]<0 ~ "SignificantDown"))
  tbb$Significant[which(!selectDiff)]="NotSignificant"
  
  tbb=as.data.frame(tbb)
  row.names(tbb)=row.names(DiffResults) #give back row.names, which is gene ID, needed for KEGGprofile
  write.csv(tbb,paste0(prefix, "_Diff.csv"))
  write.csv(tbb[selectDiff,],paste0(prefix, "_Diff_sig.csv"))
  
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

##################################
# save all comprasions designs
##################################
colnames(comprasionsTable)=c("Comparison","Factors in Comparison","Contrast")
write.csv(comprasionsTable,paste0(outFile, ".ComparisonsTable.csv"))


########################################################
# Export figures based on all results
########################################################
#vocalo plot
changeColours<-c(NotSignificant="grey",SignificantDown="blue",SignificantUp="red")
#for volcano plot
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

library(data.table)
dataForPlot=rbindlist(DiffResultsOutTableAll,idcol="Comparison")
if (useRawPvalue==1) {
  p=ggplot(dataForPlot,aes(x=logFC,y=P.Value,colour=Significant))+scale_y_continuous(trans=reverselog_trans(10),name=bquote(p~value))
} else {
  p=ggplot(dataForPlot,aes(x=logFC,y=adj.P.Val,colour=Significant))+scale_y_continuous(trans=reverselog_trans(10),name=bquote(Adjusted~p~value))
}
p=p+geom_point()+facet_wrap(~Comparison,scales="free")+scale_color_manual(values=changeColours,guide = FALSE)
saveFigure(paste0(outFile,".VolcanoPlot.png"))
print(p)
dev.off()

#venn plot of differential genes
#library("VennDiagram")
if (length(DiffResultsOutDiffGenesAll)>=1 & length(DiffResultsOutDiffGenesAll)<=5) {
  saveFigure(paste0(outFile,".Venn.png"))
  venn.diagram1(DiffResultsOutDiffGenesAll)
  dev.off()
}



save.image(paste0(outFile,".limma.RData"))














