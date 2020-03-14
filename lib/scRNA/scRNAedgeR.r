
library(edgeR)
library(ggplot2)
library(ggpubr)
library(Seurat)

sumcount<-function(ct_count, names, sample_df){
  result<-lapply(names, function(x){
    res<-ct_count[,sample_df$Cell[sample_df$Sample ==x],drop=F]
    apply(res, 1, sum)
  })
  rescount<-do.call(cbind, result)
  colnames(rescount)<-names
  return(rescount)
}

cts_cluster<-read.csv(parFile1)
cts_folder<-dirname(parFile1)

cts_unique<-rev(unique(cts_cluster$DE))

sampleGroups<-read.table(parSampleFile1, stringsAsFactors = F)
colnames(sampleGroups)<-c("Sample","Group")

comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
colnames(comparisons)<-c("Group", "Comparison")
comparisonNames<-unique(comparisons$Comparison)

ct<-cts_unique[1]
cts_name<-paste0(outFile, ".")

result<-NULL
for(ct in cts_unique){
  ct_file_name<-paste0(cts_folder, "/", cts_name, ct, ".count")
  rds_file = paste0(ct_file_name,".rds")
  sample_file = paste0(ct_file_name,".sample.csv")
  sample_df<-read.csv(sample_file, stringsAsFactors = F)
  colnames(sample_df)<-c("Cell","Sample")
  
  de_obj<-readRDS(rds_file)
  ct_count<-as.matrix(de_obj[["RNA"]]@counts)

  comp <-comparisonNames[1]
  for (comp in comparisonNames){
    prefix<-paste0(cts_name, ct, ".", comp , ".edgeR")
    comp_groups<-comparisons[comparisons$Comparison==comp,]
    
    controlGroup<-comp_groups$Group[1]
    sampleGroup<-comp_groups$Group[2]
    
    control_names<-sampleGroups$Sample[sampleGroups$Group==controlGroup]
    sample_names<-sampleGroups$Sample[sampleGroups$Group==sampleGroup]
    
    allsamples<-unique(sample_df$Sample)
    
    control_names<-control_names[control_names %in% allsamples]
    sample_names<-sample_names[sample_names %in% allsamples]
    
    if(DE_by_cell){
      cell_control<-ct_count[,sample_df$Cell[sample_df$Sample %in% control_names]]
      cell_control_group<-rep("control", ncol(cell_control))
      
      cell_sample<-ct_count[,sample_df$Cell[sample_df$Sample %in% sample_names]]
      cell_sample_group<-rep("sample", ncol(cell_sample))
    }else{
      cell_control<-sumcount(ct_count, control_names, sample_df)
      cell_control_group<-rep("control", ncol(cell_control))
      
      cell_sample<-sumcount(ct_count, sample_names, sample_df)
      cell_sample_group<-rep("sample", ncol(cell_sample))
    }
    
    cells<-cbind(cell_control, cell_sample)
    groups<-c(cell_control_group, cell_sample_group)

    dge_filename <-paste0(prefix, ".csv")
    
    #filter genes with zero count
    cells<-cells[rowSums(cells)>0,]
    
    #filter genes by tpm
    tpm = sweep(cells, 2, colSums(cells)/1e6, "/")
    min_sample<-filter_cellPercentage * ncol(cells)
    keep_rows <- rowSums(tpm > filter_minTPM) >= min_sample
    
    cells<-cells[keep_rows,]
    tpm<-tpm[keep_rows,]

    designdata<-data.frame("Group"=groups, "Cell"=colnames(cells))
    if(bComparingCluster){
      designdata$Sample=samples
    }
    
    write.csv(designdata, file=paste0(prefix, ".design"), row.names=F, quote=F)
    
    cat(prefix, "\n")
    
    dge<-DGEList(cells, group=groups)
    cat("  calcNormFactors", "\n")
    dge<-calcNormFactors(dge)
    
    cdr <- scale(colMeans(cells > 0))
    if(bComparingCluster){
      design <- model.matrix(~ cdr + samples + groups)
    }else{
      design <- model.matrix(~ cdr + groups)
    }
    rownames(design)<-colnames(cells)
    write.csv(design, file=paste0(prefix, ".design_matrix.csv"), quote=F)
    
    cat("  estimateDisp", "\n")
    dge<-estimateDisp(dge,design=design)
    
    cat("  glmQLFit", "\n")
    fitqlf<-glmQLFit(dge,design=design,robust=TRUE)
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
    
    siggenes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
    sigGenenameFile<-paste0(prefix, "_sig_genename.txt")
    write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)

    if (nrow(siggenes) > 0){
      cell_obj=de_obj[,colnames(cells)]
      geneexps=FetchData(cell_obj,vars=siggenes$gene)
      geneexps=t(geneexps)
      geneexps<-geneexps[siggenes$gene,colnames(cells)]
      sigFile<-paste0(prefix, "_sig_exp.csv")
      write.csv(cbind(sigout, geneexps), file=sigFile, quote=F)

      coords<-data.frame(cell_obj@reductions$umap@cell.embeddings)
      xlim<-c(min(coords$UMAP_1-0.1), max(coords$UMAP_1+0.1))
      ylim<-c(min(coords$UMAP_2-0.1), max(coords$UMAP_2+0.1))
      
      ddata<-designdata
      rownames(ddata)<-ddata$Cell
      cell_obj$Group=ifelse(ddata[colnames(cells), "Group"]=="control", controlGroup, sampleGroup)
      pdf(file=paste0(prefix, ".sig_genename.pdf"), onefile = T, width=21, height=7)
      siggene<-siggenes$gene[1]
      for (siggene in siggenes$gene){
        logFC<-sigout[siggene, "logFC"]
        FDR<-sigout[siggene,"FDR"]
        
        geneexp=FetchData(cell_obj,vars=c(siggene))
        colorRange<-c(min(geneexp), max(geneexp))
        fix.sc <- scale_color_gradientn(colors=c("lightgrey", "blue"), limits = colorRange)
        
        title<-paste0(siggene, ' : logFC = ', round(logFC, 2), ", FDR = ", formatC(FDR, format = "e", digits = 2))
        p0<-VlnPlot(cell_obj, feature=siggene, group.by ="Group") + NoLegend()
        if(bComparingCluster){
          p1<-DimPlot(cell_obj, reduction = "umap", label=T, group.by="seurat_clusters") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5)) + xlim(xlim) + ylim(ylim)
          p2<-FeaturePlot(object = cell_obj, features=as.character(siggene), order=T)
        }else{
          subcells<-colnames(cell_obj)[cell_obj$Group == controlGroup]
          subobj<-subset(cell_obj, cells=subcells)
          p1<-FeaturePlot(object = subobj, features=siggene, order=T) + ggtitle(paste0("Control: ", controlGroup))
          p1<-suppressMessages(expr = p1 + xlim(xlim) + ylim(ylim) + fix.sc)
          
          subcells<-colnames(cell_obj)[cell_obj$Group == sampleGroup]
          subobj<-subset(cell_obj, cells=subcells)
          p2<-FeaturePlot(object = subobj, features=siggene, order=T) + ggtitle(paste0("Sample: ", sampleGroup))
          p2<-suppressMessages(expr = p2  + xlim(xlim) + ylim(ylim) + fix.sc)
          
        }
        p<-ggarrange(plotlist=list(p0,p1,p2),nrow =1)
        g<-ggpubr::annotate_figure(
          p = p,
          top = ggpubr::text_grob(label = title, face = 'bold', size=20)
        )
        print(g)
        #break
      }
      dev.off()
    }
    
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
