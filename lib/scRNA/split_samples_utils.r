library(reshape2)
library(ggplot2)
library(Seurat)
library(gridExtra)
library(ggExtra)
library(patchwork)

is_one<-function(value, defaultValue=FALSE){
  if(is.null(value)){
    return(defaultValue)
  }
  if(is.na(value)){
    return(defaultValue)
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


rplot<-function(object, features, assay, identName, withAllCells=FALSE, n_row=1){
  DefaultAssay(object = object) <- assay
  data <- FetchData(object = object, vars = c(features, identName))
  mdata<-reshape2::melt(data, id.vars=identName)
  if (withAllCells) {
    mdata2<-mdata
    mdata2[,1] = "All cells"
    mdata<-rbind(mdata, mdata2)
  }

  nIdent=length(unique(mdata[,identName]))
  nFeature=length(features)

  gfinal=lapply(features, function(feature){
    ddata=mdata[mdata$variable==feature,]
    mvalue=ceiling(max(ddata$value))
    breaks = seq(0, mvalue, 0.5)
    
    g<-ggplot(ddata, aes(x=value)) + 
      geom_histogram(aes(y=after_stat(density)), bins=50, colour="black", fill="white", position="identity") + 
      geom_density(color="red") +
      ggtitle(feature) + 
      theme_classic()+
      theme(plot.title=element_text(hjust = 0.5),
            axis.title=element_blank(),
            axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      scale_x_continuous(breaks=breaks)
    if(nIdent > 1){
      g<-g+facet_grid(reformulate(".", identName), scale="free_y") + 
        theme(strip.background=element_blank())
    }
    return(g)
  })
  gg = wrap_plots(grobs=gfinal, nrow=n_row)

  return(gg)
}

output_tag_dist<-function(obj, filename, tagnames=rownames(obj[["HTO"]]), identName="orig.ident"){
  n_col=ceiling(sqrt(length(tagnames)))
  n_row=ceiling(length(tagnames) / n_col)
  width=min(20000, n_col * 3000 + 200)
  height=min(20000, n_row * 1500)

  g=rplot(obj, assay = "HTO", features = tagnames, identName=identName,n_row=n_row)
  ggsave(filename, g, width=width, height=height, dpi=300, units="px", bg="white")
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
  
  hto_names=unique(as.character(obj$HTO_classification))
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
  
  nFeature=length(tagnames)
  nIdent=length(unique(obj$HTO_classification))

  width=max(1800, length(tagnames) * 1500 + 300)
  
  Idents(obj) <- "HTO_classification"
  g = RidgePlot(obj, assay = "HTO", features = tagnames, ncol = length(tagnames)) &
    theme(axis.title.y=element_blank())
  ggsave(paste0(output_prefix, ".class.ridge.png"), g, width=nFeature * 4, height=(nIdent + 2) * 0.6, units="in", dpi=300, bg="white")
  
  old_levels=levels(obj$HTO_classification)
  obj$HTO_classification<-factor(obj$HTO_classification, levels=rev(old_levels))

  g = rplot(obj, assay = "HTO", features = tagnames, identName="HTO_classification")
  ggsave(paste0(output_prefix, ".class.dist.png"), g, width=nFeature * 3, height=nIdent * 1.4, units="in", dpi=300, bg="white")

  obj$HTO_classification<-factor(obj$HTO_classification, levels=old_levels)
  
  if (length(tagnames) == 2) {
    g = FeatureScatter(object = obj, feature1 = tagnames[1], feature2 = tagnames[2],group.by="HTO_classification", plot.cor=FALSE) +
      theme(aspect.ratio=1)
    ggsave(paste0(output_prefix, ".class.point.png"), g, width=5, height=3.2, units="in", dpi=300, bg="white")
  }
  
  counts=GetAssayData(obj, assay = "HTO", layer = "counts")

  tmat=data.frame(t(data.frame(counts, check.names = F)), check.names = F)
  rownames(tmat)=colnames(obj)
  tmat$HTO = unlist(obj$HTO_classification)
  tmat$HTO.global = unlist(obj$HTO_classification.global)
  write.csv(tmat, file=paste0(output_prefix, ".csv"))
  
  if(length(tagnames) > 2) {
    if(is_seurat_5_plus(obj)){
      VariableFeatures(obj) = tagnames
      cat("ScaleData...\n")
      obj<-ScaleData(obj,  features = VariableFeatures(obj))
    }else{
      obj <- FindVariableFeatures(obj, selection.method = "mean.var.plot")

      cat("ScaleData...\n")
      obj<-ScaleData(obj,  features = VariableFeatures(obj))
    }
    cat("RunUMAP...\n")
    #https://jlmelville.github.io/uwot/abparams.html
    #adjust param for umap
    obj<-RunUMAP(obj, features=tagnames, slot="scale.data", min.dist=umap_min_dist, n.neighbors=umap_num_neighbors)

    saveRDS(obj, file=paste0(output_prefix, ".umap.rds"))

    if(is_seurat_5_plus(obj)){
      umap<-FetchData(obj, c("umap_1", "umap_2")) |>
        dplyr::rename(UMAP_1 = umap_1, UMAP_2 = umap_2)
    }else{
      umap<-FetchData(obj, c("UMAP_1", "UMAP_2"))
    }
    scaled_data<-FetchData(obj, tagnames)
    colnames(scaled_data)<-paste0("Scaled_", tagnames)
    
    alldata<-cbind(umap, scaled_data, tmat)
    write.csv(alldata, file=paste0(output_prefix, ".csv"))
    
    g<-MyDimPlot(obj, reduction = "umap", group.by="HTO_classification")
    ggsave(paste0(output_prefix, ".umap.class.png"), g, width=2000, height=1800, units="px", dpi=300, bg="white")

    ncol=ceiling(sqrt(length(tagnames)))
    nrow=ceiling(length(tagnames)/ncol)
    png(paste0(output_prefix, ".umap.tag.png"), width=ncol*1500, height=1500*nrow, res=300)
    g<-FeaturePlot(obj, features=tagnames, reduction = "umap", ncol=ncol)
    ggsave(paste0(output_prefix, ".umap.tag.png"), g, width=ncol*1500, height=1500*nrow, units="px", dpi=300, bg="white")

    glist=list()
    for(tagname in tagnames){
      cols=rep("gray", length(hto_names))
      names(cols)=hto_names
      cols[tagname]="red"
      glist[[tagname]] = MyDimPlot(obj, reduction = "umap", label=T, group.by="HTO_classification", order=c(tagname))+
        scale_color_manual(values=cols) + ggtitle(tagname)
    }
    g<-wrap_plots(glist, nrow = nrow, ncol = ncol)
    ggsave(paste0(output_prefix, ".umap.assigned.png"), g, width=ncol*1500, height=1500*nrow, units="px", dpi=300, bg="white")

    cols=rep("gray", length(hto_names))
    names(cols)=hto_names
    cols[['Negative']]="blue"
    cols[["Doublet"]]="red"

    ggsave(paste0(output_prefix, ".umap.all.png"), g, width=2000, height=1800, units="px", dpi=300, bg="white")
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

do_scDemultiplex<-function(fname, rdsfile, output_prefix, init_by="demuxmix", cutoff_startval=0.5, cur_tags=NULL, cutoff_tbl=NULL ){
  library(scDemultiplex)
  library(tictoc)

  cat(fname, ":", rdsfile, " ...\n")

  obj=read_hto(rdsfile, output_prefix, cur_tags)

  p.cut=0.001

  refine_rds<-paste0(fname, ".scDemultiplex.refine.rds")

  if(init_by == "demuxmix"){
    library('demuxmix')
    library('dplyr')
    print(paste0("starting ", fname, " by demuxmix ..."))
    tic()

    hto_counts <- as.matrix(GetAssayData(obj[["HTO"]], slot = "counts"))
    dmm <- demuxmix(hto_counts, model = "naive")
    dmm_calls <- dmmClassify(dmm)
    toc1=toc()

    obj$demuxmix <- case_when(
      dmm_calls$Type == "multiplet" ~ "Doublet",
      dmm_calls$Type %in% c("negative", "uncertain") ~ "Negative",
      .default = dmm_calls$HTO)

    obj$demuxmix.global <- case_when(
      dmm_calls$Type == "multiplet" ~ "Doublet",
      dmm_calls$Type %in% c("negative", "uncertain") ~ "Negative",
      .default = "Singlet")

    obj$HTO_classification<-obj$demuxmix
    obj$HTO_classification.global<-obj$demuxmix.global

    init_column = "demuxmix";
  }else if(init_by == "HTODemux"){
    print(paste0("starting ", fname, " by HTODemux ..."))
    tic()
    obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
    toc1=toc()
    obj$HTO_classification[obj$HTO_classification.global == "Doublet"] = "Doublet"
    obj$HTODemux = obj$HTO_classification
    obj$HTODemux.global = obj$HTO_classification.global

    init_column = "HTODemux"
  }else{
    print(paste0("starting ", fname, " by cutoff ..."))
    cutoff_list = NULL
    if(!is.null(cutoff_tbl)){
      if(fname %in% cutoff_tbl$V3){
        cur_cutoff_tbl = cutoff_tbl[cutoff_tbl$V3==fname,]
        cutoff_list = split(cur_cutoff_tbl$V1, cur_cutoff_tbl$V2)
      }
    }
    tic()
    obj<-demulti_cutoff(obj, output_prefix, cutoff_startval, mc.cores=nrow(obj), cutoff_list=cutoff_list)
    toc1=toc()

    obj$HTO_classification<-obj$scDemultiplex_cutoff
    obj$HTO_classification.global<-obj$scDemultiplex_cutoff.global

    obj$scDemultiplex<-obj$scDemultiplex_cutoff
    obj$scDemultiplex.global<-obj$scDemultiplex_cutoff.global

    init_column = "scDemultiplex_cutoff";
  }

  if (nrow(obj) <= 2){
    cat(paste0("There are only ", nrow(obj), " hashtags, refinement ignored.\n"))
    toc2=NULL
  }else{
    cat(paste0("output ", init_by, " result...\n"))
    output_post_classification(obj, paste0(output_prefix, ".", init_by))

    tic(paste0("refining ", fname, " ...\n"))
    obj<-demulti_refine(obj, paste0(output_prefix, ".", init_by), init_column=init_column, mc.cores=nrow(obj), iterations = 3)
    toc2=toc()

    obj$HTO_classification<-obj$scDemultiplex
    obj$HTO_classification.global<-obj$scDemultiplex.global
  }
  saveRDS(obj, refine_rds)

  saveRDS(list("cutoff"=toc1, "refine"=toc2), paste0(output_prefix, ".scDemultiplex.tictoc.rds"))

  cat(paste0("output result...\n"))
  output_post_classification(obj, output_prefix)

  return(obj)
}
