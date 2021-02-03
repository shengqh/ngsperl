readFilesAndFormat=function(compDeseq2File) {
  deseq2<-try(read.csv(compDeseq2File, header=T, row.names=1, stringsAsFactors = F),silent=TRUE)
  if (class(deseq2)=="try-error") {
    deseq2<-try(read.csv(compDeseq2File, header=T, stringsAsFactors = F))
  }
  
  if (all(c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") %in% colnames(deseq2))) { #DESeq2 result format
    deseq2<-deseq2[,c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") ]
  } else {
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
  }
  return(result)
}


getGeneCol=function(deseq2) {
  result=list()
  if ("Feature_gene_name" %in% colnames(deseq2)) {
    result[["colName"]]="Feature_gene_name"
  } else if ("Gene" %in% colnames(deseq2)) {
    result[["colName"]]="Gene"
  }
  return(result)
}

