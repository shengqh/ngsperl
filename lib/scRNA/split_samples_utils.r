library(reshape2)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(ggExtra)

rplot<-function(object, features, assay, identName, withAllCells=FALSE){
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object, vars = c(features, identName))
  mdata<-melt(data, id.vars=identName)
  if (withAllCells) {
    mdata2<-mdata
    mdata2[,1] = "All cells"
    mdata<-rbind(mdata, mdata2)
  }
  gfinal=list()
  for(feature in features){
    ddata=mdata[mdata$variable==feature,]
    g<-ggplot(ddata, aes_string(x="value")) + 
      geom_histogram(aes(y=..density..), bins=50, colour="black", fill="white", position="identity") + 
      geom_density(color="red") +
      facet_grid(reformulate(".", identName), scale="free_y") + 
      xlab(feature) + theme_bw() + theme(strip.background=element_rect(colour="black", fill=NA),
                                         strip.text = element_text(size = 24),
                                         axis.text=element_text(size=18),
                                         axis.title=element_text(size=24))
    if (feature != features[1]){  
      g = g + ylab("")
    }
    gfinal = append(gfinal, list(g))
  }
  grid.arrange(grobs=gfinal, nrow=1)
}

read_hto<-function(h5file, output_prefix, hashtag_regex=NA) {
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
    cat("After hash tag regex filter: ", paste(rownames(htos), collapse=","), "\n")
  }else{
    htos<-mat
  }

  rowsum<-apply(htos>0, 1, sum)
  htos<-htos[rowsum > (ncol(htos) / 2),]

  cat("After zero count filter: ", paste(rownames(htos), collapse=","), "\n")

  rownames(htos)<-gsub("^TotalSeqC_", "", rownames(htos))
  rownames(htos)<-gsub("^TotalSeq_", "", rownames(htos))
  rownames(htos)<-gsub('.TotalSeqC$', "", rownames(htos))

  cat("After name clean: ", paste(rownames(htos), collapse=","), "\n")

  empty_cell_sum<-apply(htos, 2, sum)
  htos<-htos[,empty_cell_sum > 0]

  write.csv(htos, file=paste0(output_prefix, ".hto.exp.csv"))
  
  obj <- CreateSeuratObject(counts = htos, assay="HTO")
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = obj) <- "HTO"
  
  #Idents(obj) <- "HTO_classification"
  tagnames=rownames(obj[["HTO"]])
  
  width=max(1600, length(tagnames) * 1000)
  height=1400
  png(paste0(output_prefix, ".tag.dist.png"), width=width, height=height, res=300)
  rplot(obj, assay="HTO", features = tagnames, identName="orig.ident")
  dev.off()
  
  if (length(tagnames) == 2) {
    png(paste0(output_prefix, ".tag.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2], cols = "black"))
    dev.off()
  }

  return(obj)
}

output_post_classification<-function(obj, output_prefix){
  tagnames=rownames(obj[["HTO"]])

  hto_names=unique(obj$HTO_classification)
  a_hto_names=hto_names[!(hto_names %in% c("Doublet","Negative"))]
  a_hto_names=a_hto_names[order(a_hto_names)]
  hto_names=c(a_hto_names, "Negative", "Doublet")
  obj$HTO_classification=factor(obj$HTO_classification, levels=hto_names)

  width=max(1600, length(tagnames) * 1000)

  Idents(obj) <- "HTO_classification"
  png(paste0(output_prefix, ".class.ridge.png"), width=width, height=max(1400, length(tagnames) * 300), res=300)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()

  png(paste0(output_prefix, ".class.dist.png"), width=width, height=max(1400, length(tagnames) * 500), res=300)
  rplot(obj, assay = "HTO", features = tagnames, identName="HTO_classification")
  dev.off()
  
  if (length(tagnames) == 2) {
    png(paste0(output_prefix, ".class.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification"))
    dev.off()
  }
  
  tmat=data.frame(t(data.frame(obj@assays$HTO@counts)))
  rownames(tmat)=colnames(obj)
  tmat$HTO = unlist(obj$HTO_classification)
  tmat$HTO.global = unlist(obj$HTO_classification.global)
  write.csv(tmat, file=paste0(output_prefix, ".csv"))
  
  if(length(tagnames) > 2) {
    VariableFeatures(obj)<-tagnames
    obj<-ScaleData(obj)
    obj<-RunUMAP(obj, features=rownames(obj))
    
    png(paste0(output_prefix, ".umap.class.png"), width=2000, height=1800, res=300)
    g<-DimPlot(obj, reduction = "umap", group.by="HTO_classification")
    print(g)
    dev.off()
    
    nwidth=ceiling(sqrt(length(tagnames)))
    nheight=ceiling(length(tagnames)/nwidth)
    png(paste0(output_prefix, ".umap.tag.png"), width=nwidth*1500, height=1500*nheight, res=300)
    g<-FeaturePlot(obj, features=tagnames, reduction = "umap")
    print(g)
    dev.off()
    
    cols=rep("gray", length(hto_names))
    names(cols)=hto_names
    cols[['Negative']]="blue"
    cols[["Doublet"]]="red"

    png(paste0(output_prefix, ".umap.all.png"), width=2000, height=1800, res=300)
    g<-DimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification", order=c("Negative", "Doublet"))+
      scale_color_manual(values=cols)
    print(g)
    dev.off()
  }
}
