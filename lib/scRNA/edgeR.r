
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

finalList<-readRDS(parFile1)
obj<-finalList$obj

clusterDf<-read.csv(parFile2, stringsAsFactors = F, row.names=1)
obj[[cluster_name]]<-clusterDf[names(obj$orig.ident), cluster_name]

sampleGroups<-read.table(parSampleFile1, stringsAsFactors = F)
colnames(sampleGroups)<-c("Sample","Group")

comparisons<-read.table(parSampleFile2, stringsAsFactors = F)
colnames(comparisons)<-c("Group", "Comparison")
comparisonNames<-unique(comparisons$Comparison)

comp <-comparisonNames[1]

designMatrix<-NULL
for (comp in comparisonNames){
  comp_groups<-comparisons[comparisons$Comparison==comp,]
  
  controlGroup<-comp_groups$Group[1]
  sampleGroup<-comp_groups$Group[2]
  
  control_names<-sampleGroups$Sample[sampleGroups$Group==controlGroup]
  sample_names<-sampleGroups$Sample[sampleGroups$Group==sampleGroup]
  
  if(bBetweenCluster){
    prefix<-paste0(outFile, ".", comp, ".edgeR")
    
    control_cells<-rownames(clusterDf)[clusterDf[,cluster_name] %in% control_names]  
    sample_cells<-rownames(clusterDf)[clusterDf[,cluster_name] %in% sample_names]  
  
    all_cells<-c(control_cells, sample_cells)

    de_obj<-subset(obj, cells=all_cells)
    de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
    de_obj$DisplayGroup<-c(rep(controlGroup, length(control_cells)), rep(sampleGroup, length(sample_cells)))
      
    designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
    designfile<-paste0(prefix, ".design")
    write.csv(designdata, file=designfile, row.names=F, quote=F)
    
    curdf<-data.frame(prefix=prefix, cellType="", comparison=comp, sampleInGroup=1, design=designfile, stringsAsFactors = F)
    if (is.null(designMatrix)){
      designMatrix = curdf
    }else{
      designMatrix = rbind(designMatrix, curdf)
    }
  }else{
    cts = unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name])
    prefixList<-gsub(" ", "_", cts)
    prefixList<-gsub(":", "_", prefixList)
    prefixList<-gsub("_+", "_", prefixList)
    
    idx<-1
    for (idx in c(1:length(cts))){
      ct = cts[idx]
      prefix = paste0(outFile, ".", prefixList[idx], ".", comp, ".edgeR")
      
      clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
      de_obj<-subset(obj, cells=clusterCt$cell)
      clusterCt$sample=de_obj$orig.ident
      
      control_cells<-rownames(clusterCt)[clusterCt$sample %in% control_names]  
      sample_cells<-rownames(clusterCt)[clusterCt$sample %in% sample_names]  
      
      all_cells<-c(control_cells, sample_cells)
      
      de_obj<-subset(obj, cells=all_cells)
      de_obj$Group<-c(rep("control", length(control_cells)), rep("sample", length(sample_cells)))
      de_obj$DisplayGroup<-c(rep(controlGroup, length(control_cells)), rep(sampleGroup, length(sample_cells)))
      
      designdata<-data.frame("Group"=de_obj$Group, "Cell"=colnames(de_obj), "Sample"=de_obj$orig.ident, "DisplayGroup"=de_obj$DisplayGroup)
      designfile<-paste0(prefix, ".design")
      write.csv(designdata, file=designfile, row.names=F, quote=F)
      
      curdf<-data.frame(prefix=prefix, cellType=ct, comparison=comp, sampleInGroup=0, design=designfile, stringsAsFactors = F)
      if (is.null(designMatrix)){
        designMatrix = curdf
      }else{
        designMatrix = rbind(designMatrix, curdf)
      }
    }
  }
}

result<-NULL
idx<-1
for(idx in c(1:nrow(designMatrix))){
  prefix=designMatrix[idx, "prefix"]
  designfile=designMatrix[idx, "design"]
  cellType=designMatrix[idx, "cellType"]
  comp=designMatrix[idx, "comparison"]
  sampleInGroup=designMatrix[idx, "sampleInGroup"]
  
  designdata<-read.csv(designfile, stringsAsFactors = F)
  groups<-designdata$Group
  
  de_obj<-subset(obj, cells=designdata$Cell)
  cells<-as.matrix(de_obj[["RNA"]]@counts)
    
  #filter genes with zero count
  cells<-cells[rowSums(cells)>0,]
  
  #filter genes by tpm
  tpm = sweep(cells, 2, colSums(cells)/1e6, "/")
  min_sample<-filter_cellPercentage * ncol(cells)
  keep_rows <- rowSums(tpm > filter_minTPM) >= min_sample
  rm(tpm)
  
  cells<-cells[keep_rows,]
  cdr <- scale(colMeans(cells > 0))
  
  if(sampleInGroup){
    samples<-designdata$Sample
    design <- model.matrix(~ cdr + samples + groups)
  }else{
    design <- model.matrix(~ cdr + groups)
  }
  
  rownames(design)<-colnames(cells)
  write.csv(design, file=paste0(prefix, ".design_matrix.csv"), quote=F)
  
  cat(prefix, "\n")
  
  dge<-DGEList(cells, group=groups)
  cat("  calcNormFactors", "\n")
  dge<-calcNormFactors(dge)

  cat("  estimateDisp", "\n")
  dge<-estimateDisp(dge,design=design)
  
  cat("  glmQLFit", "\n")
  fitqlf<-glmQLFit(dge,design=design,robust=TRUE)
  qlf<-glmQLFTest(fitqlf)
  out<-topTags(qlf, n=Inf)
  dge_filename <-paste0(prefix, ".csv")
  write.csv(out$table, file=dge_filename, quote=F)

  if(useRawPvalue){
    sigout<-out$table[(out$table$PValue<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }else{
    sigout<-out$table[(out$table$FDR<=pvalue) & (abs(out$table$logFC)>=log2(foldChange)),]
  }
  sigFile<-paste0(prefix, ".sig.csv")
  write.csv(sigout, file=sigFile, quote=F)
  
  siggenes<-data.frame(gene=rownames(sigout), stringsAsFactors = F)
  sigGenenameFile<-paste0(prefix, ".sig_genename.txt")
  write.table(siggenes, file=sigGenenameFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  gseaFile<-paste0(prefix, "_GSEA.rnk")
  rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * out$table$F)
  write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)
  
  curDF<-data.frame("prefix"=prefix, "cellType"=cellType, "comparison"=comp, "betweenCluster"=bBetweenCluster, "sampleInGroup"=sampleInGroup, "deFile"=dge_filename, "sigFile"=sigFile, "sigGenenameFile"=sigGenenameFile, "gseaFile"=gseaFile, "designFile"=designfile)
  if(is.null(result)){
    result<-curDF
  }else{
    result<-rbind(result, curDF)
  }
}

write.csv(result, file=paste0(outFile, ".edgeR.files.csv"), quote=F)
