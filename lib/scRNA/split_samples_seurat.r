library(Seurat)
library(ggplot2)

split<-function(h5file, output_prefix, hashtag_regex=NA) {
  if(grepl('.rds$', h5file)){
    meta<-readRDS(h5file)
    meta<-as.matrix(meta)
    #tag number should less than cell number
    if(ncol(meta) < nrow(meta)){
      meta=t(meta)
    }
  }else{
    sdata<-Read10X_h5(h5file)
    meta<-sdata[[2]]
    meta<-as.matrix(meta)
  }
  mat<-meta
  write.csv(mat, file=paste0(output_prefix, ".alltags.exp.csv"))

  cat("All tags: ", paste(rownames(mat), collapse=","), "\n")
  
  if (!is.na(hashtag_regex)) {
    htos<-mat[grepl(hashtag_regex, rownames(mat)),]
    if (nrow(htos) == 0){
      stop(paste0("Cannot find hashtag based on regex ", hashtag_regex, " for tags ", paste(rownames(mat), collapse=",")))
    }
    cat("All hash tags matched regex: ", paste(rownames(htos), collapse=","), "\n")
  }else{
    htos<-mat
  }
  rownames(htos)<-gsub("^TotalSeqC_", "", rownames(htos))
  rownames(htos)<-gsub("^TotalSeq_", "", rownames(htos))
  rownames(htos)<-gsub('.TotalSeqC$', "", rownames(htos))

  empty_cell_sum<-apply(htos, 2, sum)
  htos<-htos[,empty_cell_sum > 0]

  write.csv(htos, file=paste0(output_prefix, ".hto.exp.csv"))

  pbmc.hashtag <- CreateSeuratObject(counts = htos, assay="HTO")
  pbmc.hashtag <- NormalizeData(pbmc.hashtag, normalization.method = "CLR")
  pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
  pbmc.hashtag$HTO_classification[pbmc.hashtag$HTO_classification.global == "Doublet"] = "Doublet"

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
  write.csv(tmat, paste0(output_prefix, ".csv"))
  
  VariableFeatures(pbmc.hashtag)<-tagnames
  pbmc.hashtag<-ScaleData(pbmc.hashtag)
  pbmc.hashtag<-RunUMAP(pbmc.hashtag, features=rownames(pbmc.hashtag))
  
  png(paste0(output_prefix, ".umap.class.png"), width=1000, height=800)
  g<-DimPlot(pbmc.hashtag, reduction = "umap", group.by="HTO_classification")
  print(g)
  dev.off()
  
  png(paste0(output_prefix, ".umap.tag.png"), width=1600, height=1600)
  g<-FeaturePlot(pbmc.hashtag, features=tagnames, reduction = "umap")
  print(g)
  dev.off()
  
  hto_names=unique(pbmc.hashtag$HTO_classification)
  a_hto_names=hto_names[!(hto_names %in% c("Doublet","Negative"))]
  a_hto_names=a_hto_names[order(a_hto_names)]
  hto_names=c(a_hto_names, "Negative", "Doublet")
  cols=rep("gray", length(hto_names))
  names(cols)=hto_names
  cols[['Negative']]="blue"
  cols[["Doublet"]]="red"

  pbmc.hashtag$HTO_classification=factor(pbmc.hashtag$HTO_classification, levels=hto_names)
  png(paste0(output_prefix, ".umap.all.png"), width=1000, height=800)
  g<-DimPlot(pbmc.hashtag, reduction = "umap", label=T, group.by="HTO_classification", order=c("Negative", "Doublet"))+
    scale_color_manual(values=cols)
  print(g)
  dev.off()

}

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  h5file = "/data/cqs/seurat_data/hto12_hto_valid.rds"
  output_prefix = "/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_cutoff/result/hto12/hto12.HTO"
  hashtag_regex='Hashtag|TotalSeqC_|C025|Benign|Tumor|HTO|HEK|THP|K562|KG1'
}else{
  h5file = args[1]
  output_prefix = args[2]
  hashtag_regex = args[3]
}

print(paste0("h5file=", h5file))
print(paste0("output_prefix=", output_prefix))
print(paste0("hashtag_regex=", hashtag_regex))

split(h5file, output_prefix, hashtag_regex)