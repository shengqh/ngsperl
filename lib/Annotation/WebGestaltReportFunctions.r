
updateEnrichedTable<-function(enriched){
	rowCount<-nrow(enriched)

	enriched$tblShortCaption<-apply(enriched, 1, function(x){tabRef(paste0("enriched ",x["geneSet"]))})
	enriched$tblLongCaption<-apply(enriched, 1, function(x){tabRef(paste0("enriched ",x["geneSet"]), paste0("Significantly differential expressed genes in ",x["description"]))})
	enriched$tblLink<-tolower(gsub(' ','', enriched$tblShortCaption))

	for(idx in c(1:rowCount)){
		entry<-enriched[idx,]
		enriched$link[idx]<-paste0("[", entry$tblShortCaption, "](#", entry$tblLink, ")")
	}

	return(enriched)
}

readFilesAndFormat=function(compDeseq2File) {
  library(data.table)
  library(dplyr)
  library(tibble)
  deseq2<-data.frame(fread(compDeseq2File, header=T), row.names=1)
  
  deseq2_columns = c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange")
  if (all(deseq2_columns %in% colnames(deseq2))) { #DESeq2 result format
    deseq2<-deseq2[,deseq2_columns ]
  }else if(all(c("logFC", "logCPM", "PValue", "FDR") %in% colnames(deseq2))) { #edgeR result format
    deseq2 = deseq2 %>% tibble::rownames_to_column(var="GENE") %>%
              dplyr::mutate(FoldChange = 2**logFC) %>%
              dplyr::select(GENE, everything())              
  }
  return(deseq2)
}

getDiffCol=function(deseq2) {
  result=list()
  if ("FoldChange" %in% colnames(deseq2)) {
    result[["colName"]]="FoldChange"
    result[["centerValue"]]=1
  } else if ("meth.diff" %in% colnames(deseq2)) {
    result[["colName"]]="meth.diff"
    result[["centerValue"]]=0
  } else {
    stop("Can't find diff column")
  }
  return(result)
}

getGeneCol=function(deseq2) { #try to find gene column
  result=list()
  
  result[["colName"]]=-1
  if (ncol(deseq2)==1) {
    result[["colName"]]=1
  } else if ("Feature_gene_name" %in% colnames(deseq2)) { #deseq2 result format
    result[["colName"]]="Feature_gene_name"
  } else if ("Hugo_Symbol" %in% colnames(deseq2)) { #maf file
    result[["colName"]]="Hugo_Symbol"
  } else { #other formats
    geneInd=grep("Gene|gene|GENE",colnames(deseq2))
    if (length(geneInd)>0) { #first column name with "gene"
      result[["colName"]]=geneInd[1]
    } else { 
      colClass<-sapply(deseq2, class)
      countNotNumIndex<-which((colClass!="numeric" & colClass!="integer"))
      if(length(countNotNumIndex)>0) {
        result[["colName"]]=countNotNumIndex[1]
      }else{
        result[["colName"]]=-1
      }
    }
  }
  return(result)
}
