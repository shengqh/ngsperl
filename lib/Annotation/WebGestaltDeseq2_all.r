rm(list=ls()) 
outFile='P9061'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_harmony_dr0.2_nrh_03_choose_edgeR_inCluster_byCell_WebGestalt/result/P9061.WebGestaltR.files.csv'
parFile2='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_harmony_dr0.2_nrh_03_choose_edgeR_inCluster_byCell/result/P9061.edgeR.files.csv'
parFile3=''
outputDirectory='.'


setwd('/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_harmony_dr0.2_nrh_03_choose_edgeR_inCluster_byCell_WebGestalt_link_edgeR/result')

### Parameter setting end ###

source("WebGestaltReportFunctions.r")
library("rmarkdown")

annoFiles<-read.csv(parFile1, header=T, stringsAsFactors = F)

deseq2Files<-read.csv(parFile2, header=T, stringsAsFactors = F, row.names=1)
deseq2Files$deFile<-file.path(dirname(parFile2), deseq2Files$deFile)

comparisons<-unique(annoFiles$Comparison)
comparison<-comparisons[1]
for (comparison in comparisons){
  compAnnoFiles<-annoFiles$File[annoFiles$Comparison == comparison]
  compDeseq2File<-deseq2Files$deFile[deseq2Files$prefix == comparison]
  
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
