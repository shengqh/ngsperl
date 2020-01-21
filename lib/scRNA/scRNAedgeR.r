
library(edgeR)
library(heatmap3)

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

drawHCA<-function(prefix, rldselect, ispaired, designData, conditionColors, gnames, outputFormat, usePearsonInHCA=TRUE){
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
      labRow=row.names(rldselect)
      margins=c(12,8)
    } else {
      labRow=NA
      margins=c(12,5)
    }
    
    filePrefix<-paste0(prefix, ".heatmap")
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
                 main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),  
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
                 main=paste0("Hierarchical Cluster Using ", genecount, " Genes"),  
                 cexCol=cexCol, 
                 useRaster=FALSE,
                 legendfun=function() showLegend(legend=paste0("Group ", gnames), col=c("red","blue"),cex=1.0,x="center"))
      }
      dev.off()
    }
  }
}

sumcount<-function(ct_count, names, sample_df){
  result<-lapply(names, function(x){
    res<-ct_count[,sample_df$Cell[sample_df$Sample ==x],drop=F]
    apply(res, 1, sum)
  })
  rescount<-do.call(cbind, result)
  colnames(rescount)<-names
  return(rescount)
}

comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
comparisonNames<-unique(comparisons$V3)

cts_cluster<-read.csv(parFile1)
cts_folder<-dirname(parFile1)
cts_unique<-rev(unique(cts_cluster$DE))

ct<-cts_unique[1]
cts_name<-paste0(outFile, ".")

result<-NULL
for(ct in cts_unique){
  ct_file_name<-paste0(cts_folder, "/", cts_name, ct, ".count")
  rds_file = paste0(ct_file_name,".rds")
  sample_file = paste0(ct_file_name,".sample.csv")
  
  ct_count<-readRDS(rds_file)
  sample_df<-read.csv(sample_file)
  colnames(sample_df)<-c("Cell","Sample")
  
  comp <-comparisonNames[1]
  for (comp in comparisonNames){
    prefix<-paste0(cts_name, ct, ".", comp , ".edgeR")
    dge_filename <-paste0(prefix, ".csv")
    
    #if(file.exists(dge_filename)){
    #  next
    #}

    comp_groups<-comparisons[comparisons$V3==comp,]
    control_names<-comp_groups$V1[comp_groups$V2=="control"]
    sample_names<-comp_groups$V1[comp_groups$V2=="sample"]

    control_names<-control_names[control_names %in% sample_df$Sample]
    sample_names<-sample_names[sample_names %in% sample_df$Sample]
    
    if(DE_by_cell){
      cell_control<-ct_count[,sample_df$Cell[sample_df$Sample %in% control_names]]
      cell_control_group<-rep(1, ncol(cell_control))
      
      cell_sample<-ct_count[,sample_df$Cell[sample_df$Sample %in% sample_names]]
      cell_sample_group<-rep(2, ncol(cell_sample))
    }else{
      cell_control<-sumcount(ct_count, control_names, sample_df)
      cell_control_group<-rep(1, ncol(cell_control))
      
      cell_sample<-sumcount(ct_count, sample_names, sample_df)
      cell_sample_group<-rep(2, ncol(cell_sample))
    }

    cells<-cbind(cell_control, cell_sample)
    groups<-c(cell_control_group, cell_sample_group)

    #filter genes with zero count
    cells<-cells[rowSums(cells)>0,]
    
    #filter genes by tpm
    tpm = sweep(cells, 2, colSums(cells)/1e6, "/")
    min_sample<-filter_cellPercentage * ncol(cells)
    keep_rows <- rowSums(tpm > filter_minTPM) >= min_sample
    
    cells<-cells[keep_rows,]
    tpm<-tpm[keep_rows,]

    designdata<-data.frame("Group"=groups, "Sample"=colnames(cells))
    write.csv(designdata, file=paste0(prefix, ".design"), row.names=F, quote=F)
    
    cat(prefix, "\n")
    
    dge<-DGEList(cells, group=groups)
    cat("  calcNormFactors", "\n")
    dge<-calcNormFactors(dge)
    
    cdr <- scale(colMeans(cells > 0))
    design <- model.matrix(~ cdr + groups)
    
    cat("  estimateDisp", "\n")
    dge<-estimateDisp(dge,design=design)
    
    cat("  glmQLFit", "\n")
    fitqlf<-glmQLFit(dge,design=design)
    qlf<-glmQLFTest(fitqlf)
    out<-topTags(qlf, n=Inf)
    outTpm<-tpm[rownames(out$table),]
    write.csv(cbind(out$table, outTpm), file=dge_filename, quote=F)

    if(useRawPvalue){
      sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
    }else{
      sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
    }
    sigTpm<-tpm[rownames(sigout),]
    sigFile<-paste0(prefix, "_sig.csv")
    write.csv(cbind(sigout, sigTpm), file=sigFile, quote=F)
    
    siggenes<-data.frame(gene=rownames(sigout))
    sigGenenameFile<-paste0(prefix, "_sig_genename.txt")
    write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
    
    gseaFile<-paste0(prefix, "_GSEA.rnk")
    rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * out$table$F)
    write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)
    
    #cat("  heatmap", "\n")
    #y <- cpm(dge, log=TRUE)
    #conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[group]))
    #gnames<-c(control_name, sample_name)
    #drawHCA(prefix, y, FALSE, NULL, conditionColors, gnames, "png", usePearsonInHCA=TRUE)
    
    curDF<-data.frame("celltype"=ct, "comparison"=comp, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile)
    if(is.null(result)){
      result<-curDF
    }else{
      result<-rbind(result, curDF)
    }
  }
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)