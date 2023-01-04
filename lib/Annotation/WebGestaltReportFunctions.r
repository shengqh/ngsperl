
readFilesAndFormat=function(compDeseq2File) {
  if (grepl(".csv$",basename(compDeseq2File))) { #csv file
    readFun<-function(...) read.csv(...)
  } else {
    readFun<-function(...) read.table(...,sep="\t")
  }
  deseq2<-try(readFun(compDeseq2File, header=T, row.names=1, stringsAsFactors = F),silent=TRUE)
  if (class(deseq2)=="try-error") {
    deseq2<-try(readFun(compDeseq2File, header=T, stringsAsFactors = F))
  }
  
  if (all(c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") %in% colnames(deseq2))) { #DESeq2 result format
    deseq2<-deseq2[,c("Feature_gene_name", "baseMean", "pvalue", "padj", "FoldChange") ]
  }else if(all(c("logFC", "logCPM", "F", "PValue", "FDR") %in% colnames(deseq2))) {
    deseq2$GENE<-rownames(deseq2)
    deseq2$FoldChange=2**deseq2$logFC
    deseq2=deseq2[c("GENE", "logFC", "logCPM", "F", "PValue", "FDR", "FoldChange")]
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


getGeneCol=function(deseq2) { #try to find gene column
  result=list()
  
  if (ncol(deseq2)==1) {
    result[["colName"]]=1
  } else if ("Feature_gene_name" %in% colnames(deseq2)) { #deseq2 result format
    result[["colName"]]="Feature_gene_name"
  } else if ("Hugo_Symbol" %in% colnames(deseq2)) { #maf file
    result[["colName"]]="Hugo_Symbol"
  } else { #otehr formats
    geneInd=grep("Gene|gene|GENE",colnames(deseq2))
    if (length(geneInd)>0) { #first column name with "gene"
      result[["colName"]]=geneInd[1]
    } else { 
      #guess gene column, by contents with both number and character (so not all numeric data, can be gene IDs)
      #also want column with many unique values (so that won't be columns like sample names)
      geneInd=which.max(apply(geneList[-1,],2,function(x) 
        min(length(intersect(grep("[a-zA-Z][a-zA-Z]",x),grep("[0-9]",x))),length(unique(x)))
        ))
      result[["colName"]]=geneInd[1]
    }
  }
  return(result)
}
