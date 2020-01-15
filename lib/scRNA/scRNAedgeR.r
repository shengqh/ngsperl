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

comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
comparisonNames<-unique(comparisons$V3)

cts_cluster<-read.csv(parFile1)
cts_folder<-dirname(parFile1)
cts_unique<-rev(unique(cts_cluster$DE))

ct<-"20_T_cells"
cts_name<-paste0(outFile, ".")

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
    
    if(file.exists(dge_filename)){
      next
    }

    comp_groups<-comparisons[comparisons$V3==comp,]
    control_names<-comp_groups$V1[comp_groups$V2=="control"]
    sample_names<-comp_groups$V1[comp_groups$V2=="sample"]
  
    cell_control<-ct_count[,sample_df$Cell[sample_df$Sample %in% control_names]]
    cell_control_group<-rep(1, ncol(cell_control))
  
    cell_sample<-ct_count[,sample_df$Cell[sample_df$Sample %in% sample_names]]
    cell_sample_group<-rep(2, ncol(cell_sample))

    cat(prefix, "\n")
    
    dge<-DGEList(cbind(cell_control, cell_sample), group = c(cell_control_group, cell_sample_group))
    cat("  calcNormFactors", "\n")
    dge<-calcNormFactors(dge)
    
    group<-dge$samples[,"group"]
    design<-model.matrix(~group)
    cat("  estimateDisp", "\n")
    dge<-estimateDisp(dge,design)
    
    cat("  glmQLFit", "\n")
    fitqlf<-glmQLFit(dge,design)
    qlf<-glmQLFTest(fitqlf,coef=2)
    out<-topTags(qlf, n=Inf, adjust.method = "BH")
    
    write.csv(out, file=dge_filename)
    
    sigout<-data.frame(gene=rownames(out)[(out$table$FDR<=0.05) & (abs(out$table$logFC)>=1)])
    write.table(sigout, file=paste0(prefix, "_sig_genename.txt"), row.names=F, col.names=F, sep="\t", quote=F)
    
    rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * out$table$F)
    write.table(rankout, file=paste0(prefix, "_GSEA.rnk"), row.names=F, col.names=F, sep="\t", quote=F)
    
    #cat("  heatmap", "\n")
    #y <- cpm(dge, log=TRUE)
    #conditionColors<-as.matrix(data.frame(Group=c("red", "blue")[group]))
    #gnames<-c(control_name, sample_name)
    #drawHCA(prefix, y, FALSE, NULL, conditionColors, gnames, "png", usePearsonInHCA=TRUE)
  }
}

