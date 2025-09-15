library(Seurat)
library(ggplot2)

split<-function(h5file, output_prefix, hashtag_regex=NA) {
  sdata<-Read10X_h5(h5file)

  exp<-sdata[[1]]
  meta<-sdata[[2]]
  mat<-as.matrix(meta)
  rowsum<-apply(mat>0, 1, sum)
  mat<-mat[rowsum > (ncol(mat) / 2),]
  rownames(mat)<-gsub("\\.1$", "", rownames(mat))
  
  if (!is.na(hashtag_regex)) {
    htos<-mat[grepl(hashtag_regex, rownames(mat)),]
  }else{
    htos<-mat
  }
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
  
  
  #Idents(pbmc.hashtag) <- "HTO_classification"
  tagnames=rownames(pbmc.hashtag[["HTO"]])
  
  Idents(pbmc.hashtag) <- "orig.ident"
  width=max(10, length(tagnames) * 5)
  pdf(paste0(output_prefix, ".tag.dist.pdf"), width=width, height=8)
  print(RidgePlot(pbmc.hashtag, assay = "HTO", features = tagnames, ncol = length(tagnames), cols = "black", fill="white"))
  dev.off()

  if (length(tagnames) == 2) {
    pdf(paste0(output_prefix, ".tag.point.pdf"), width=width, height=8)
    print(FeatureScatter(object = pbmc.hashtag, feature1 = tagnames[1], feature2 = tagnames[2], cols = "black"))
    dev.off()
  }

  Idents(pbmc.hashtag) <- "HTO_classification"
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
  g<-ggplot(tmat, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + facet_grid(rows="HTO.global") + theme_bw() + theme(strip.background = element_blank())
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

if (length(args) == 0) {
  h5file = "C:/Users/sheng/projects/paula_hurley/20201208_scRNA_split/filtered_feature_bc_matrix.h5"
  output_prefix = "C:/Users/sheng/projects/paula_hurley/20201208_scRNA_split/split_samples/result/HYW_4701/HYW_4701.HTO"
  hashtag_regex = NA
}else{
  h5file = args[1]
  output_prefix = args[2]
  hashtag_regex = args[3]
}

print(paste0("h5file=", h5file))
print(paste0("output_prefix=", output_prefix))
print(paste0("hashtag_regex=", hashtag_regex))

split(h5file, output_prefix, hashtag_regex)