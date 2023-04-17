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

  #empty_cell_sum<-apply(htos, 2, sum)
  #htos<-htos[,empty_cell_sum > 0]

  write.csv(htos, file=paste0(output_prefix, ".hto.exp.csv"))

  obj <- CreateSeuratObject(counts = htos, assay="HTO")
  obj <- NormalizeData(obj, assay="HTO", normalization.method = "CLR")
  obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
  
  #RidgePlot(obj, assay = "HTO", features = c("HEK-A", "K562-B", "KG1-A", "THP1-C"), ncol = 2)
  
  obj$HTO_classification[obj$HTO_classification.global == "Doublet"] = "Doublet"

  #Idents(obj) <- "HTO_classification"
  tagnames=rownames(obj[["HTO"]])
  
  Idents(obj) <- "orig.ident"
  width=max(10, length(tagnames) * 5)
  pdf(paste0(output_prefix, ".tag.dist.pdf"), width=width, height=8)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames), cols = "black", fill="white"))
  dev.off()

  if (length(tagnames) == 2) {
    pdf(paste0(output_prefix, ".tag.point.pdf"), width=width, height=8)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2], cols = "black"))
    dev.off()
  }

  Idents(obj) <- "HTO_classification"
  width=max(10, length(tagnames) * 5)
  pdf(paste0(output_prefix, ".dist.pdf"), width=width, height=8)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()
  
  tmat=data.frame(t(htos))
  tmat$HTO = obj$HTO_classification
  tmat$HTO.global = obj$HTO_classification.global
  write.csv(tmat, paste0(output_prefix, ".csv"))
  
  VariableFeatures(obj)<-tagnames
  obj<-ScaleData(obj)
  obj<-RunUMAP(obj, features=rownames(obj))
  
  png(paste0(output_prefix, ".umap.class.png"), width=1000, height=800)
  g<-MyDimPlot(obj, reduction = "umap", group.by="HTO_classification")
  print(g)
  dev.off()
  
  png(paste0(output_prefix, ".umap.tag.png"), width=1600, height=1600)
  g<-FeaturePlot(obj, features=tagnames, reduction = "umap")
  print(g)
  dev.off()
  
  hto_names=unique(obj$HTO_classification)
  a_hto_names=hto_names[!(hto_names %in% c("Doublet","Negative"))]
  a_hto_names=a_hto_names[order(a_hto_names)]
  hto_names=c(a_hto_names, "Negative", "Doublet")
  cols=rep("gray", length(hto_names))
  names(cols)=hto_names
  cols[['Negative']]="blue"
  cols[["Doublet"]]="red"

  obj$HTO_classification=factor(obj$HTO_classification, levels=hto_names)
  png(paste0(output_prefix, ".umap.all.png"), width=1000, height=800)
  g<-MyDimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification", order=c("Negative", "Doublet"))+
    scale_color_manual(values=cols)
  print(g)
  dev.off()

}

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  h5file = "/data/cqs/seurat_data/hto12_hto_valid.rds"
  output_prefix = "/scratch/cqs/shengq2/papers/20210703_scrna_hto/hto_samples_HTODemux/result/hto12/hto12.HTO"
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