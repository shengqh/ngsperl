library(reshape2)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(ggExtra)

rplot<-function(object, features, assay, identName, withAllCells=FALSE, n_row=1){
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
      xlab(feature) + theme_bw()+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=24),
            axis.title.y=element_blank())
    if(length(unique(mdata[,identName])) > 1){
      g<-g+facet_grid(reformulate(".", identName), scale="free_y") + 
        theme(strip.background=element_rect(colour="black", fill=NA),
              strip.text = element_text(size = 24))
    }
    gfinal = append(gfinal, list(g))
  }
  grid.arrange(grobs=gfinal, nrow=n_row)
}

read_hto<-function(rdsfile, output_prefix) {
  htos<-readRDS(rdsfile)

  obj <- CreateSeuratObject(counts = htos, assay="HTO")
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
  DefaultAssay(object = obj) <- "HTO"

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
    obj<-ScaleData(obj, assay="HTO")
    obj<-RunUMAP(obj, assay="HTO", slot="scale.data", features=rownames(obj))
    
    saveRDS(object = obj, file = paste0(output_prefix, ".umap.rds"))
    
    png(paste0(output_prefix, ".umap.class.png"), width=2000, height=1800, res=300)
    g<-DimPlot(obj, reduction = "umap", group.by="HTO_classification")
    print(g)
    dev.off()
    
    nwidth=ceiling(sqrt(length(tagnames)))
    nheight=ceiling(length(tagnames)/nwidth)
    png(paste0(output_prefix, ".umap.tag.png"), width=nwidth*1500, height=1500*nheight, res=300)
    g<-FeaturePlot(obj, features=tagnames, reduction = "umap", ncol=nwidth)
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

build_summary<-function(allres, output_prefix, nameMapFile=NA){
  colnames(allres)<-c("File", "Sample")
  dat=lapply(allres$File, function(x){
    dd=read.csv(x, check.names=F)
    table(dd$HTO.global)
  })
  dat.all=do.call(rbind, dat)
  rownames(dat.all)=allres$Sample
  write.csv(dat.all, file=paste0(output_prefix, ".summary.csv"), quote=F)
  
  mdat=reshape2::melt(dat.all)
  colnames(mdat)=c("Sample", "Class", "Cell")
  
  png(paste0(output_prefix, ".summary.png"), width=1600, height=1200, res=300)
  g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(g)
  dev.off()
}
