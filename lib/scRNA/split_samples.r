library(Seurat)
library(ggplot2)

split<-function(h5file, output_prefix) {
  sdata<-Read10X_h5(h5file)

  exp<-sdata[[1]]
  meta<-sdata[[2]]
  mat<-as.matrix(meta)
  rowsum<-apply(mat>0, 1, sum)
  mat<-mat[rowsum > (ncol(mat) / 2),]
  rownames(mat)<-gsub("\\.1$", "", rownames(mat))
  
  htos<-mat[grepl("Hash", rownames(mat)),]
  rownames(htos)<-gsub("_.*", "", rownames(htos))
  
  pbmc.hashtag <- CreateSeuratObject(counts = exp)
  pbmc.hashtag <- NormalizeData(pbmc.hashtag)
  # Find and scale variable features
  pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
  pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
  
  pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = htos)
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
  pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
  
  table(pbmc.hashtag$HTO_classification.global)
  
  Idents(pbmc.hashtag) <- "HTO_maxID"
  tagnames=rownames(pbmc.hashtag[["HTO"]])
  
  width=max(10, length(tagnames) * 5)
  pdf(paste0(output_prefix, ".dist.pdf"), width=width, height=8)
  print(RidgePlot(pbmc.hashtag, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()
  
  tmat=data.frame(t(mat))
  tmat$HTO = pbmc.hashtag$HTO_classification
  tmat$HTO.global = pbmc.hashtag$HTO_classification.global
  
  tmat=data.frame(t(mat))
  tmat$HTO = pbmc.hashtag$HTO_classification
  tmat$HTO.global = pbmc.hashtag$HTO_classification.global
  tmat$nCount_RNA = pbmc.hashtag$nCount_RNA
  tmat$nFeature_RNA = pbmc.hashtag$nFeature_RNA
  write.csv(tmat, paste0(output_prefix, ".csv"))
  
  width=length(tagnames) * 3
  pdf(paste0(output_prefix, ".class.pdf"), width=width, height=4)
  g<-ggplot(tmat, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + facet_grid(~HTO.global) + theme_bw() + theme(strip.background = element_blank())
  print(g)
  dev.off()
  
  pbmc.hashtag$hash.ID = factor(pbmc.hashtag$hash.ID, levels=sort(unique(as.character(pbmc.hashtag$hash.ID))))
  Idents(pbmc.hashtag) <- "hash.ID"

  pdf(paste0(output_prefix, ".nCount.pdf"), width=10, height=8)
  print(VlnPlot(pbmc.hashtag, features="nCount_RNA"))
  dev.off()
  
  pdf(paste0(output_prefix, ".nFeature.pdf"), width=10, height=8)
  print(VlnPlot(pbmc.hashtag, features="nFeature_RNA"))
  dev.off()
}

args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
outputPrefix = args[2]

split(inputFile, outputPrefix)