###############################################################
#make heatmap for genes in each category of WebGestalt result
###############################################################

usePearsonInHCA<-TRUE
hmcols <- colorRampPalette(c("green", "black", "red"))(256)
maxCategory=20 #max 20 heatmap (categories)
if (!exists("fileTypePrefix")) { #to add in file name, GOID+fileTypePrefix.png
  fileTypePrefix="_DESeq2-vsd-heatmap" 
}

library(heatmap3)

#function from DESeq2.r
drawHCA<-function(prefix, rldselect, ispaired, designData, conditionColors, gnames, outputFormat,
                  main=paste0("Hierarchical Cluster Using ", nrow(rldselect), " Genes"),
                  labRow=ifelse(nrow(rldselect)<=30,row.names(rldselect),NA),
                  fileTypePrefix="_DESeq2-vsd-heatmap"){
  genecount<-nrow(rldselect)
  showRowDendro = genecount <= 50
  if(genecount > 2){
    cexCol = max(1.0, 0.2 + 1/log10(ncol(rldselect)))
    if(ispaired){
      htColors<-rainbow(length(unique(designData$Paired)))
      gsColors<-as.matrix(data.frame(Group=conditionColors, Sample=htColors[designData$Paired]))
    }else{
      gsColors = conditionColors;
    }
    if (genecount<=30) {
      #labRow=row.names(rldselect)
      margins=c(12,8)
    } else {
      #labRow=NA
      margins=c(12,5)
    }
    
    filePrefix<-paste0(prefix, fileTypePrefix)
    for(format in outputFormat){
      openPlot(filePrefix, format, 10, 10, 3000, 3000, "HCA")
      if(usePearsonInHCA){
        heatmap3(rldselect, 
                 col = hmcols, 
                 ColSideColors = gsColors, 
                 margins=margins, 
                 scale="r", 
                 labRow=labRow,
                 showRowDendro=showRowDendro,
                 main=main,  
                 cexCol=cexCol, 
                 useRaster=FALSE,
                 legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
      }else{
        heatmap3(rldselect, 
                 col = hmcols, 
                 ColSideColors = gsColors, 
                 margins=margins, 
                 scale="r", 
                 distfun=dist, 
                 labRow=labRow,
                 showRowDendro=showRowDendro,
                 main=main,  
                 cexCol=cexCol, 
                 useRaster=FALSE,
                 legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
      }
      dev.off()
    }
  }
}

openPlot<-function(filePrefix, format, pdfWidth, pdfHeight, otherWidth, otherHeight, figureName){
  fileName<-paste0(filePrefix, ".", tolower(format))
  if(format == "PDF"){
    pdf(fileName, width=pdfWidth, height=pdfHeight, useDingbats=FALSE)
  }else if(format == "TIFF"){
    tiff(filename=fileName, width=otherWidth, height=otherHeight, res=300)
  }else {
    png(filename=fileName, width=otherWidth, height=otherHeight, res=300)
  }
  cat("saving", figureName, "to ", fileName, "\n")
}



annoFiles<-read.table(parSampleFile1, header=F, sep="\t", stringsAsFactors = F)
deseq2Files<-read.table(parSampleFile2, header=F, sep="\t", stringsAsFactors = F)
deseq2VsdFiles<-read.table(parSampleFile3, header=F, sep="\t", stringsAsFactors = F)
deseq2DesignFiles<-read.table(parSampleFile4, header=F, sep="\t", stringsAsFactors = F)

comparisons<-unique(annoFiles$V2)
comparison<-comparisons[1]
for (comparison in comparisons){
  dir.create(comparison,showWarnings=FALSE)
  compAnnoFiles<-annoFiles[annoFiles$V2 == comparison,]$V1
  compDeseq2File<-deseq2Files[deseq2Files$V2 == comparison,1]
  deseq2VsdFile<-deseq2VsdFiles[deseq2VsdFiles$V2 == comparison,1]
  deseq2DesignFile<-deseq2DesignFiles[deseq2DesignFiles$V2 == comparison,1]
  
  compAnnoFile<-compAnnoFiles[1]
  for(compAnnoFile in compAnnoFiles){
    if (!file.exists(compAnnoFile)) {
      warning(paste0(compAnnoFile," doesn't exist, Skip!"))
      next;
    }
    targetFolder=paste0("./",comparison,"/",basename(compAnnoFile),"/")
    dir.create(targetFolder,showWarnings=FALSE)
    
    currentWd=getwd()
    setwd(targetFolder)
    
    category <- gsub(paste0(".*?", comparison,"_"), "", basename(compAnnoFile) )
    category <- gsub(".txt", "", category )
    category <- gsub("_", " ", category )
    
    enriched<-read.table(compAnnoFile, sep="\t", header=T, stringsAsFactors = F)
    #deseq2<-read.csv(compDeseq2File, header=T, row.names=1, stringsAsFactors = F)
    #deseq2<-deseq2[,c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") ]
    deseq2=readFilesAndFormat(compDeseq2File)
    geneCol=getGeneCol(deseq2)[["colName"]]
    if(geneCol == -1){
      samples=colnames(deseq2)
      deseq2$gene_name=rownames(deseq2)
      deseq2=deseq2[,c("gene_name", samples)]
      geneCol="gene_name"
    }
    diffCol=getDiffCol(deseq2)[["colName"]]
    diffCenterValue=getDiffCol(deseq2)[["centerValue"]]
    
    #geneExpression<-read.csv(deseq2VsdFile, header=T, row.names=1, stringsAsFactors = F)
    geneExpression<-readFilesAndFormat(deseq2VsdFile)
    
    designData<-read.delim(deseq2DesignFile, header=T,  stringsAsFactors = F,check.names=F)
    conditionColors=as.matrix(data.frame(Group=c("red", "blue")[as.factor(designData[,2])]))
    gnames=levels(as.factor(designData[,2]))
    
    rowCount<-min(maxCategory,nrow(enriched))
    idx<-1
    for(idx in c(1:rowCount)){
      entry<-enriched[idx,]
      userIds<-unlist(strsplit( entry$userId[1], ';'))
      #entryTable<-geneExpression[row.names(deseq2)[which(deseq2$Feature_gene_name %in% userIds)],]
      #row.names(entryTable)=make.unique(deseq2[which(deseq2$Feature_gene_name %in% userIds),"Feature_gene_name"])
      entryTable<-geneExpression[row.names(deseq2)[which(deseq2[,geneCol] %in% userIds)],]
      row.names(entryTable)=make.unique(deseq2[which(deseq2[,geneCol] %in% userIds),geneCol])
      
      drawHCA(prefix=make.names(enriched[idx,1]), rldselect=entryTable, ispaired=FALSE, designData=deseq2Design, 
                   conditionColors=conditionColors, gnames=gnames, outputFormat="PNG",
              main=paste(enriched[idx,1:2],collapse=", "),labRow=row.names(entryTable),fileTypePrefix=fileTypePrefix)
        
    }
    setwd(currentWd)
  }
}


