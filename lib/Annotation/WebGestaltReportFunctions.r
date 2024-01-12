
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
