library("rmarkdown")

annoFiles<-read.csv(parFile1, stringsAsFactors = F)

deseq2Files<-read.csv(parFile2, stringsAsFactors = F)
deseq2Files$deFile=paste0(dirname(parFile2), "/", deseq2Files$deFile)

comparisons<-unique(annoFiles$Comparison)
comparison<-comparisons[1]
for (comparison in comparisons){
  subFiles=subset(annoFiles, annoFiles$Comparison == comparison)
  compAnnoFiles<-annoFiles$File[annoFiles$Comparison == comparison]
  compDeseq2File<-deseq2Files$deFile[deseq2Files$comparison == comparison]
  
  for(sidx in c(1:nrow(subFiles))) {
    compAnnoFile=subFiles$File[sidx]
    category=subFiles$Database[sidx]

    if (!file.exists(compAnnoFile)) {
      warning(paste0(compAnnoFile," doesn't exist, Skip!"))
      next;
    }

    category <- gsub("_", " ", category )
    
    enriched<-read.table(compAnnoFile, sep="\t", header=T, stringsAsFactors = F)
    deseq2=readFilesAndFormat(compDeseq2File)
    #deseq2<-read.csv(compDeseq2File, header=T, row.names=1, stringsAsFactors = F)
    #deseq2<-deseq2[,c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") ]
    diffCol=getDiffCol(deseq2)[["colName"]]
    diffCenterValue=getDiffCol(deseq2)[["centerValue"]]
    geneCol=getGeneCol(deseq2)[["colName"]]
    
    rowCount<-nrow(enriched)
    idx<-1
    for(idx in c(1:rowCount)){
	    entry<-enriched[idx,]
	    userIds<-unlist(strsplit( entry$userId[1], ';'))
	    #entryTable<-deseq2[deseq2$Feature_gene_name %in% userIds,]
	    #geneUp<-sum(entryTable$FoldChange > 1)
	    #geneDown<-sum(entryTable$FoldChange < 1)
	    entryTable<-deseq2[deseq2[,geneCol] %in% userIds,]
	    if (!is.null(diffCol)) {
	      geneUp<-sum(entryTable[,diffCol] > diffCenterValue)
	      geneDown<-sum(entryTable[,diffCol] < diffCenterValue)
	      enriched$geneUp[idx]<-geneUp
	      enriched$geneDown[idx]<-geneDown
	    }
  		enriched$geneSet[idx]<-paste0("[", entry$geneSet, "](", entry$link, ")")
  		enriched$link[idx]<-""
    }

    plotData <- list(enriched = enriched,
                     deseq2 = deseq2,
                     category = category)
    
    output_path <- paste0(normalizePath("."), "/", basename(compAnnoFile), ".html")
    saveRDS(plotData, paste0(output_path, ".rds"))
    
    #output_dir = "E:/temp"
    #output_file = "temp1.html"
    output_dir = dirname(output_path)
    output_file = basename(output_path)
    
    cat("Output report to:", output_path, "\n")
    rmarkdown::render("WebGestaltDeseq2.rmd",
                      output_dir = output_dir,
                      output_file = output_file,
                      params = list(data = plotData))
  }
}
