rm(list=ls()) 
outFile='P11389_SLAMseq'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''
outputDirectory='.'


setwd('/nobackup/brown_lab/projects/20240520_11389_SLAMseq_mm10/deseq2_tc_read_sizeFactor_WebGestalt_link_deseq2/result')

### Parameter setting end ###

source("WebGestaltReportFunctions.r")
library("rmarkdown")

if (!exists("annoFiles")) {
  annoFiles<-read.table(parSampleFile1, header=F, sep="\t", stringsAsFactors = F)
}
if (!exists("deseq2Files")) {
  deseq2Files<-read.table(parSampleFile2, header=F, sep="\t", stringsAsFactors = F)
}

comparisons<-unique(annoFiles$V2)
comparison<-comparisons[1]
for (comparison in comparisons){
  compAnnoFiles<-annoFiles[annoFiles$V2 == comparison,]$V1
  compDeseq2File<-deseq2Files[deseq2Files$V2 == comparison,1]
  
  compAnnoFile<-compAnnoFiles[1]
  for(compAnnoFile in compAnnoFiles){
    if (!file.exists(compAnnoFile)) {
      warning(paste0(compAnnoFile," doesn't exist, Skip!"))
      next;
    }
    category <- gsub(paste0(".*?", comparison,"_"), "", basename(compAnnoFile) )
    category <- gsub(".txt", "", category )
    category <- gsub("_", " ", category )
    
    enriched<-read.table(compAnnoFile, sep="\t", header=T, stringsAsFactors = F)
    deseq2=readFilesAndFormat(compDeseq2File)
    #deseq2<-read.csv(compDeseq2File, header=T, row.names=1, stringsAsFactors = F)
    #deseq2<-deseq2[,c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") ]
    geneCol=getGeneCol(deseq2)[["colName"]]
    if(geneCol == -1){
      samples=colnames(deseq2)
      deseq2$gene_name=rownames(deseq2)
      deseq2=deseq2[,c("gene_name", samples)]
      geneCol="gene_name"
    }
    diffCol=getDiffCol(deseq2)[["colName"]]
    diffCenterValue=getDiffCol(deseq2)[["centerValue"]]
    
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
  		enriched$geneSet[idx]<-paste0("[", entry$geneSet, "](", entry$link, "){target='_blank'}")
  		enriched$link[idx]<-""
    }

    plotData <- list(enriched = enriched,
                     deseq2 = deseq2,
                     category = category)
    
    fname<-gsub("^enrichment_results_", "",  basename(compAnnoFile))
    fname<-gsub('.txt$', '', fname)
    output_path <- paste0(normalizePath("."), "/", fname, ".html")
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

