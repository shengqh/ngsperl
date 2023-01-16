library(reshape2)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(ggExtra)

is_one<-function(value){
  if(is.na(value)){
    return(FALSE)
  }
  if(is.null(value)){
    return(FALSE)
  }
  return(value == '1')
}

read_map<-function(map_file){
  tbl=read.table(map_file, sep="\t")
  result=split(tbl$V1, tbl$V2)
}

save_to_matrix<-function(counts, target_folder) {
  if(!dir.exists(target_folder)){
    dir.create(target_folder)
  }
  
  bar_file=paste0(target_folder, "/barcodes.tsv")
  writeLines(colnames(counts), bar_file)
  gzip(bar_file, overwrite=T)
  
  feature_file=paste0(target_folder, "/features.tsv")
  writeLines(rownames(counts), feature_file)
  gzip(feature_file, overwrite=T)
  
  matrix_file=paste0(target_folder, "/matrix.mtx")
  writeMM(counts, matrix_file)
  gzip(matrix_file, overwrite=T)
}

output_tag_dist<-function(obj, filename){
  tagnames=rownames(obj[["HTO"]])
  n_col=ceiling(sqrt(length(tagnames)))
  n_row=ceiling(length(tagnames) / n_col)
  width=min(20000, n_col * 3000 + 200)
  height=min(20000, n_row * 1500)
  png(filename, width=width, height=height, res=300)
  rplot(object=obj, assay="HTO", features = tagnames, identName="orig.ident", n_row=n_row)
  dev.off()
}

# rename_tags<-function(tags){
#   result<-gsub("^TotalSeqC_", "", tags)
#   result<-gsub("^TotalSeq_", "", result)
#   result<-gsub('.TotalSeqC$', "", result)
#   return(result)
# }

#read raw hto file
read_hto_file<-function(cname, cfiles){
  if(length(cfiles) == 1){
    cat("preparing", cname, ":", cfiles, " ...\n")
    if(grepl('.h5$', cfiles)){
      sdata<-Read10X_h5(cfiles)
      exp<-sdata[[1]]
      writeLines(colnames(exp),con = paste0(cname, ".barcodes.tsv"))
      htos<-sdata[[2]]
      htos<-as.matrix(htos)
    }else{
      stop(paste0("Unknown file format:", cfiles))
    }
  }else{
    hfile=cfiles[[1]]
    cat("preparing", cname, ":", hfile, " ...\n")
    efile=cfiles[[2]]
    if(grepl('hto_mtx.rds$', hfile)){
      htos<-readRDS(hfile)
      htos<-as.matrix(htos)
      #tag number should less than cell number
      if(ncol(htos) < nrow(htos)){
        htos=t(htos)
      }
      exp<-readRDS(efile)
      cells.use <- intersect(colnames(exp), colnames(htos))
      exp<-exp[, cells.use]
      htos<-htos[, cells.use]
    }else{
      stop(paste0("Unknown file format:", cfiles))
    }
  }

  cat("Original tag names: ", paste(rownames(htos), collapse=","), "\n")

  # cat("Before name clean: ", paste(rownames(htos), collapse=","), "\n")
  # rownames(htos)<-rename_tags(rownames(htos))
  # cat("After name clean: ", paste(rownames(htos), collapse=","), "\n")

  return(list(exp=exp, htos=htos))
}

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
    mvalue=ceiling(max(ddata$value))
    breaks = seq(0, mvalue, 0.5)
    
    g<-ggplot(ddata, aes_string(x="value")) + 
      geom_histogram(aes(y=..density..), bins=50, colour="black", fill="white", position="identity") + 
      geom_density(color="red") +
      xlab(feature) + theme_bw()+
      theme(axis.text=element_text(size=18),
            axis.title=element_text(size=24),
            axis.title.y=element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_x_continuous(breaks=breaks)
    if(length(unique(mdata[,identName])) > 1){
      g<-g+facet_grid(reformulate(".", identName), scale="free_y") + 
        theme(strip.background=element_rect(colour="black", fill=NA),
              strip.text = element_text(size = 24))
    }
    gfinal = append(gfinal, list(g))
  }
  grid.arrange(grobs=gfinal, nrow=n_row)
}

#read prepared hto object
read_hto<-function(rdsfile, output_prefix, cur_tags=NULL) {
  htos<-readRDS(rdsfile)
  if(!all(is.null(cur_tags))){
    if(!all(cur_tags %in% rownames(htos))){
      stop(paste0("Not all tags ", paste0(cur_tags, collapse = ","), " in tag names: ", paste0(rownames(htos), collapse = ",")))
    }
  }

  if('Seurat' %in% class(htos)){
    obj<-htos
    if(!all(is.null(cur_tags))){
      obj<-subset(obj, features=cur_tags)
    }
  }else{
    if(!all(is.null(cur_tags))){
      htos<-htos[cur_tags,]
    }
    obj <- CreateSeuratObject(counts = htos, assay="HTO")
    obj <- subset(obj, nCount_HTO > 0)
    # Normalize HTO data, here we use centered log-ratio (CLR) transformation
    obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
    DefaultAssay(object = obj) <- "HTO"
  }

  return(obj)
}

output_post_classification<-function(obj, output_prefix, umap_min_dist=0.3, umap_num_neighbors=30, tagnames=NULL){
  if(all(is.null(tagnames))){
    tagnames=rownames(obj[["HTO"]])
  }
  
  hto_names=unique(obj$HTO_classification)
  a_hto_names=hto_names[!(hto_names %in% c("Doublet","Negative"))]
  a_hto_names=a_hto_names[order(a_hto_names)]
  hto_names=c(a_hto_names, "Negative", "Doublet")
  obj$HTO_classification=factor(obj$HTO_classification, levels=hto_names)

  tbl<-as.data.frame(table(obj$HTO_classification))
  colnames(tbl)<-c("Tagname", "Cell")

  if(length(unique(obj$orig.ident)) > 1){
    sdata<-unique(FetchData(obj, c("orig.ident", "HTO_classification")))
    if(nrow(sdata) == length(hto_names)){
      name_map = unlist(split(sdata$orig.ident, sdata$HTO_classification))
      tbl$Sample = unlist(name_map[tbl$Tagname])
    }
  }
  cat_file=paste0(output_prefix, ".tag_cell.csv")
  write.csv(tbl, cat_file)
  
  width=max(1800, length(tagnames) * 1500 + 300)
  
  Idents(obj) <- "HTO_classification"
  png(paste0(output_prefix, ".class.ridge.png"), width=width, height=max(1400, (length(tagnames) + 2) * 300), res=300)
  print(RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)))
  dev.off()
  
  png(paste0(output_prefix, ".class.dist.png"), width=width, height=max(1400, (length(tagnames) + 2) * 500), res=300)
  old_levels=levels(obj$HTO_classification)
  obj$HTO_classification<-factor(obj$HTO_classification, levels=rev(old_levels))
  rplot(obj, assay = "HTO", features = tagnames, identName="HTO_classification")
  obj$HTO_classification<-factor(obj$HTO_classification, levels=old_levels)
  dev.off()
  
  if (length(tagnames) == 2) {
    png(paste0(output_prefix, ".class.point.png"), width=2000, height=1800, res=300)
    print(FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification"))
    dev.off()
  }
  
  tmat=data.frame(t(data.frame(obj@assays$HTO@counts, check.names = F)), check.names = F)
  rownames(tmat)=colnames(obj)
  tmat$HTO = unlist(obj$HTO_classification)
  tmat$HTO.global = unlist(obj$HTO_classification.global)
  write.csv(tmat, file=paste0(output_prefix, ".csv"))
  
  if(length(tagnames) >= 2) {
    VariableFeatures(obj)<-tagnames

    cat("ScaleData...\n")
    obj<-ScaleData(obj)

    cat("RunUMAP...\n")
    #https://jlmelville.github.io/uwot/abparams.html
    #adjust param for umap
    obj<-RunUMAP(obj, features=tagnames, slot="scale.data", min.dist=umap_min_dist, n.neighbors=umap_num_neighbors)

    saveRDS(obj, file=paste0(output_prefix, ".umap.rds"))

    umap<-FetchData(obj, c("UMAP_1", "UMAP_2"))
    scaled_data<-FetchData(obj, tagnames)
    colnames(scaled_data)<-paste0("Scaled_", tagnames)
    
    alldata<-cbind(umap, scaled_data, tmat)
    write.csv(alldata, file=paste0(output_prefix, ".csv"))
    
    png(paste0(output_prefix, ".umap.class.png"), width=2000, height=1800, res=300)
    g<-DimPlot(obj, reduction = "umap", group.by="HTO_classification")
    print(g)
    dev.off()

    ncol=ceiling(sqrt(length(tagnames)))
    nrow=ceiling(length(tagnames)/ncol)
    png(paste0(output_prefix, ".umap.tag.png"), width=ncol*1500, height=1500*nrow, res=300)
    g<-FeaturePlot(obj, features=tagnames, reduction = "umap", ncol=ncol)
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

lookup_cutoff<-function(cutoffs, fname, tagname){
  if(is.null(cutoffs)){
    return(0)
  }
  
  fc<-cutoffs[cutoffs$V3==fname & cutoffs$V2==tagname,]
  if(nrow(fc) == 0){
    return(0)
  }
  
  return(fc$V1[1])
}

draw_cutoff<-function(prefix, values, cut_off){
  png(paste0(prefix, ".cutoff.png"), width=2000, height=1600, res=300)
  his<-hist(values,200,F,xlab="concentration",ylab="density", main=NULL,col="grey")
  lines(density(values),lwd=1.5,col="blue")
  abline(v=cut_off,lwd=1.5,col="red")

  minb=min(his$breaks)
  maxb=max(his$breaks)
  x=minb + (maxb-minb) * 3 /4
  
  y=max(his$density)
  legend(x=x, y=y, legend=c("density", "cutoff"), col=c("blue", "red"), lty=c(1,2))
  dev.off()
}

