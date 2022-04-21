library(harmony)
library(cowplot)

file_not_empty<-function(filename){
  if(is.null(filename)){
    return(FALSE)
  }
  
  if(is.na(filename)){
    return(FALSE)
  }
  
  if(!file.exists(filename)){
    return(FALSE)
  }
  
  return(file.info(filename)$size != 0)
}

calc_weight<-function(cellType){
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  return(weight)
}

read_cell_markers_file<-function(panglao5_file, species, remove_subtype_of="", HLA_panglao5_file="", curated_markers_file=""){
  #preparing cell activity database
  marker<-data.frame(fread(panglao5_file))
  if(remove_subtype_of != ""){
    remove_subtype_of<-unique(unlist(strsplit(remove_subtype_of, ",")))
    pangdb_ct <- read.table(HLA_panglao5_file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    remove_subtype_of<-remove_subtype_of[remove_subtype_of %in% marker$cell.type]
    layer="Layer2"
    removed<-c()
    for(layer in c("Layer1","Layer2", "Layer3")){
      layer_ct<-unique(pangdb_ct[,layer])
      for(rs in remove_subtype_of){
        if(rs %in% removed){
          next
        }
        
        if(rs %in% layer_ct){
          subdb<-rownames(pangdb_ct)[pangdb_ct[,layer]==rs]
          if(rs %in% subdb){ 
            subdb<-subdb[!(subdb %in% remove_subtype_of)]
            marker<-marker[!(marker$cell.type %in% subdb),]
          }else{
            marker$cell.type[marker$cell.type %in% subdb] = rs
          }
          removed<-c(removed, rs)
        }
      }
    }
  }
  
  hsind<-regexpr(species,marker$species)
  marker_species<-marker[hsind>0 & marker$ubiquitousness.index<0.05,]
  if (species=="Mm") {
    ##change the gene symbol only keep the first letter capitalize
    marker_species$official.gene.symbol<-toMouseGeneSymbol(marker_species$official.gene.symbol)
  }
  if (species=="Hs") {
    marker_species$official.gene.symbol<-toupper(marker_species$official.gene.symbol)
  }
  cellType<-tapply(marker_species$official.gene.symbol,marker_species$cell.type,list)
  weight<-calc_weight(cellType)
  
  if(!is.null(curated_markers_file) && curated_markers_file != ""){
    curated_markers_df<-read.table(curated_markers_file, sep="\t", header=F, stringsAsFactors=F)
    curated_markers_celltype<-split(curated_markers_df$V2, curated_markers_df$V1)
    for(cmct in names(curated_markers_celltype)){
      cellType[[cmct]]=curated_markers_celltype[[cmct]]
    }
    weight=calc_weight(cellType)
  } 
  return(list(cellType=cellType, weight=weight))
}

get_selfdefinedCelltype <- function(file, finalLayer="layer3"){
  file.ori <- scan(file,what = "",sep="\n")
  marker.ori <- matrix("undefined",nrow=length(file.ori)/2,ncol = 2)
  for (i in 1:length(file.ori)){if(i%%2==1){marker.ori[(i+1)/2,1]=substr(file.ori[i],2,nchar(file.ori[i]))}else{marker.ori[i/2,2]=file.ori[i]}}
  ref <- as.data.frame(marker.ori[,1])
  ref <- tidyr::separate(ref,col=colnames(ref),into = c("celltype","subtypeOf"),sep = ",")
  rownames(ref) <- ref$celltype
  celltype.ori <- sapply(tapply(marker.ori[,2],marker.ori[,1],list), function(x) strsplit(x,","))
  
  celltype.tem <- as.data.frame(matrix(unlist(celltype.ori)))
  celltype.tem[,2] <- rep(names(celltype.ori),lengths(celltype.ori))
  celltype.tem <- tidyr::separate(celltype.tem,col=V2,into = c("celltype","subtypeOf"),sep = ",")
  
  ct <- as.data.frame(matrix("Undefined",nrow = nrow(celltype.tem),ncol = length(unique(celltype.tem$subtypeOf))+2),stringsAsFactors=FALSE)
  ct[,1:2] <- celltype.tem[,1:2]
  for (i in 3:ncol(ct)){
    for (j in 1:nrow(ct)) {if(ct[j,i-1] %in% rownames(ref)){ct[j,i]<-ref[ct[j,i-1],2]}else{ct[j,i]<-ct[j,i-1]}}
  }
  
  layer <- ncol(ct)
  if(! identical(ct[,ncol(ct)],ct[,ncol(ct)-1])){layer <- ncol(ct)}else{
    for (i in 1:ncol(ct)) {if(identical(ct[,i],ct[,i+1])){layer <- i;break}}
  }
  
  celltype.trim <- ct[,1:layer]
  for (i in 1:nrow(celltype.trim)){
    tag <- which(celltype.trim[i,]==celltype.trim[i,ncol(celltype.trim)])
    if(length(tag)>1){
      if(ncol(celltype.trim)+1-length(tag)==3){suf <- ct[i,3]}
      else{suf <- ct[i,3:(ncol(celltype.trim)+1-length(tag))]}
      pre <- rep(celltype.trim[i,2],length(tag)-1)
      celltype.trim[i,] <- c(celltype.trim[i,][1:2],pre,suf)
    }
  }
  
  for (i in 1:ncol(celltype.trim)){if(i == 1){colnames(celltype.trim)[i]="gene"}else{colnames(celltype.trim)[i]=paste0("layer",(ncol(celltype.trim)+1-i),"")}}
  
  ct.final <- celltype.trim
  for (i in 1:ncol(ct.final)){
    if(i >1){ct.final[,i]=celltype.trim[,ncol(celltype.trim)+2-i];colnames(ct.final)[i] <- colnames(celltype.trim)[ncol(celltype.trim)+2-i]}
  }
  
  cellType<-tapply(ct.final$gene,ct.final[[finalLayer]],list)
  weight<-calc_weight(cellType)
  return(list(cellType=cellType, weight=weight))
}

toMouseGeneSymbol<-function(x){
  result=paste0(toupper(substr(x,1,1)),tolower(substr(x,2,nchar(x))))
  return(result)
}

read_cell_cluster_file<-function(fileName, sort_cluster_name="seurat_clusters", display_cluster_name=sort_cluster_name, prefix=""){
  result<-read.csv(fileName, stringsAsFactors = F, row.names = 1)
  
  if(!display_cluster_name %in% colnames(result)){
    result[,display_cluster_name] = paste0(prefix, result[,sort_cluster_name])
  }
  
  if(sort_cluster_name=="seurat_clusters"){
    clusters = result[!duplicated(result$seurat_clusters),]
    clusters = clusters[order(clusters$seurat_clusters),]
  }else{
    clusters = result[!duplicated(result[,c("seurat_clusters", sort_cluster_name)]),]
    clusters = clusters[order(clusters[,sort_cluster_name], clusters$seurat_clusters),]
  }
  
  cluster_names=colnames(clusters)[grepl("_clusters", colnames(clusters))]
  
  for(cluster_name in cluster_names){
    cf_levels=unique(clusters[,cluster_name])
    result[,cluster_name] = factor(result[,cluster_name], levels=cf_levels)
  }
  return(result)
}

find_markers<-function(object, by_sctransform, ident.1, ident.2, min.pct = 0.5, logfc.threshold = 0.6){
  assay=ifelse(by_sctransform, "SCT", "RNA")
  markers=FindMarkers(object, assay=assay, ident.1=ident.1, ident.2=ident.2, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  markers=markers[markers$p_val_adj < 0.01,]
  return(markers)
}

find_all_markers<-function(object, by_sctransform, min.pct = 0.5, logfc.threshold = 0.6){
  assay=ifelse(by_sctransform, "SCT", "RNA")
  markers=FindAllMarkers(object, assay=assay, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  markers=markers[markers$p_val_adj < 0.01,]
  return(markers)
}

get_cluster_count<-function(counts, clusters){
  if(is.null(levels(clusters))){
    allClusters=unique(clusters)
  }else{
    allClusters=levels(clusters)
  }
  cluster<-allClusters[1]
  csums=lapply(allClusters, function(cluster){
    #cat(cluster, "\n")
    cells=names(clusters)[clusters==cluster]
    if(length(cells) == 0){
      stop("no cells")
    }
    subcounts=counts[,cells,drop=F]
    #cat(ncol(subcounts), "\n")
    Matrix::rowSums(subcounts)
  })
  
  result=do.call(cbind, csums)
  colnames(result)=allClusters
  
  gcount=Matrix::rowSums(result)
  result=result[gcount > 0,]
  return(result)
}

get_group_count=function(curobj, groupName="active.ident") {
  counts=GetAssayData(curobj,assay="RNA",slot="counts")
  curgroups=curobj[[groupName]]
  clusters=curgroups[,1]
  names(clusters)=rownames(curgroups)
  result=get_cluster_count(counts, clusters)
  return(result)
}

add_cluster<-function(object, new.cluster.name, new.cluster.ids){
  seurat_clusters<-object[["seurat_clusters"]]$seurat_clusters
  names(new.cluster.ids) <- levels(seurat_clusters)
  
  new.cluster.values<-plyr::mapvalues(x = seurat_clusters, from = levels(seurat_clusters), to = new.cluster.ids)
  names(new.cluster.values)<-names(seurat_clusters)
  
  object[[new.cluster.name]]<-new.cluster.values
  object
}

add_celltype<-function(obj, celltype_df, celltype_column){
  new.cluster.ids<-split(celltype_df[,celltype_column], celltype$seurat_clusters)
  obj[[celltype_column]] = unlist(new.cluster.ids[as.character(obj$seurat_clusters)])
  
  celltype_df$seurat_column = paste0(celltype_df$seurat_clusters, " : ", celltype_df[,celltype_column])
  new.cluster.ids<-split(celltype_df$seurat_column, celltype$seurat_clusters)
  obj[[paste0("seurat_", celltype_column)]] = unlist(new.cluster.ids[as.character(obj$seurat_clusters)])
  return(obj)
}

run_pca<-function(object, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform){
  if (by_sctransform) {
    object <- RunPCA(object = object, verbose=FALSE)
  }else{
    if (Remove_Mt_rRNA) {
      rRNA.genes <- grep(pattern = rRNApattern,  rownames(object), value = TRUE)
      Mt.genes<- grep (pattern= Mtpattern,rownames(object), value=TRUE )
      var.genes <- dplyr::setdiff(VariableFeatures(object), c(rRNA.genes,Mt.genes))
    } else {
      var.genes <- VariableFeatures(object)
    }
    object <- RunPCA(object = object, features = var.genes, verbose=FALSE)
  }
  return(object)
}

run_cluster_only<-function(object, pca_dims, resolution, random.seed, reduction="pca"){
  object <- FindNeighbors(object = object, reduction=reduction, dims=pca_dims, verbose=FALSE)
  object <- FindClusters(object=object, verbose=FALSE, random.seed=random.seed, resolution=resolution)
  return(object)
}

run_cluster<-function(object, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, min.pct = 0.5, logfc.threshold = 0.6, reduction="pca"){
  object=run_cluster_only(object, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, min.pct, logfc.threshold, reduction)
  if (by_sctransform) {
    markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  }else{
    markers <- FindAllMarkers(object, features=var.genes, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  }
  markers <- markers[markers$p_val_adj < 0.01,]
  return(list(object=object, markers=markers))
}

ORA_celltype<-function(medianexp,cellType,weight){
  ORA_result<-matrix(NA, nrow=length(cellType),ncol=dim(medianexp)[2])
  CTA_result<-matrix(0,nrow=length(cellType),ncol=dim(medianexp)[2])
  
  colnames(ORA_result)=colnames(medianexp)
  colnames(CTA_result)=colnames(medianexp)
  
  exp_z<-scale(medianexp)
  genenames<-rownames(medianexp)   
  for (j in 1: dim(medianexp)[2]){
    clusterexp<-medianexp[,j] 
    clusterexp_z<-exp_z[,j]
    for (i in 1:length(cellType)){
      
      ct_exp<-length(intersect(genenames[clusterexp>0],cellType[[i]]))
      ct_not_exp<-length(cellType[[i]])-ct_exp
      exp_not_ct<-sum(clusterexp>0)-ct_exp
      not_exp_not_ct<-length(clusterexp)-ct_not_exp 
      cont.table<-matrix(c(ct_exp,ct_not_exp,exp_not_ct,not_exp_not_ct),nrow=2)
      ORA_result[i,j]<-fisher.test(cont.table,alternative="greater")$p.value
      ###
      weight_ss<-weight[names(weight)%in%cellType[[i]]]
      ind<-match(names(weight_ss),genenames)
      exp_ss<-clusterexp_z[ind[!is.na(ind)]]
      weight_ss<-weight_ss[!is.na(ind)]
      CTA_result[i,j]<-sum(exp_ss*weight_ss)/(length(exp_ss)^(1/3))
    }
  }
  rownames(ORA_result)<-rownames(CTA_result)<-names(cellType)
  minp_ora_ind<- apply(ORA_result,2,function(x){which.min(x)})
  minp_ora<-apply(ORA_result,2,min)
  names(minp_ora)<-rownames(ORA_result)[minp_ora_ind]
  
  max_cta_ind<- apply(CTA_result,2,function(x){which.max(x)})
  max_cta<-apply(CTA_result,2,max,na.rm=T)
  names(max_cta)<-rownames(CTA_result)[max_cta_ind]
  return(list(ora=ORA_result,cta=CTA_result,min_ora=minp_ora,max_cta=max_cta))
}

get_cta_ora_mat<-function(predict_celltype) {
  cta_index<-apply(predict_celltype$cta,2,function(x){return(order(x,decreasing=T)[1:2])})
  cta_index<-unique(sort(cta_index))
  
  cta_mat<- predict_celltype$cta[cta_index,]
  colnames(cta_mat)<-paste0(names(predict_celltype$max_cta), " : ", colnames(predict_celltype$cta))
  
  ora_mat<- predict_celltype$ora[cta_index,]
  ora_mat<--log10(ora_mat)  
  colnames(ora_mat)<-colnames(cta_mat)
  
  return(list(cta_mat=cta_mat, ora_mat=ora_mat))
}

get_cta_combined<-function(obj, predicted){
  cta_table<-data.frame(Cluster=colnames(predicted$cta), 
                        CellType=names(predicted$max_cta),
                        CtaScore=round(predicted$max_cta * 10) / 10.0,
                        stringsAsFactors = F)
  cta_table$OraPvalue=apply(cta_table, 1, function(x){
    ct=x[2]
    cl=x[1]
    predicted$ora[ct, cl]
  })
  
  cluster_sample<-as.data.frame.matrix(table(obj$seurat_clusters, obj$orig.ident))
  cluster_sample<-cluster_sample[as.character(cta_table$Cluster),]
  
  nc<-apply(cluster_sample, 2, function(x){
    tc=sum(x)
    perc<-x/tc
    return(round(perc*1000) / 10.0)
  })
  colnames(nc)<-paste0(colnames(nc), "_perc")
  
  cta_combined<-cbind(cta_table, cluster_sample, nc)
  
  return(cta_combined)
}

get_batch_samples<-function(batch_file, all_sample_names){
  if(file.exists(batch_file)){
    pools = read.table(batch_file, header=F, stringsAsFactors = F)
    missnames<-all_sample_names[!(all_sample_names %in% pools$V1)]
    if(length(missnames) > 0){
      pools<-rbind(pools, data.frame(V1=missnames, V2=missnames))
    }
    poolmap=split(pools$V2, pools$V1)
  }else{
    poolmap=split(all_sample_names, all_sample_names)
  }
  return(poolmap)
}

draw_dimplot<-function(mt, filename, split.by) {
  nSplit = length(unique(mt[,split.by]))
  nWidth=ceiling(sqrt(nSplit))
  nHeight=ceiling(nSplit / nWidth)
  
  png(filename, width=nWidth * 600 + 200, height=nHeight*600, res=300)
  g1<-ggplot(mt, aes(x=UMAP_1,y=UMAP_2)) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    facet_wrap(as.formula(paste("~", split.by)), ncol=nWidth, nrow=nHeight) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  print(g1)
  dev.off()
  return(g1)
}

do_harmony<-function(obj, npcs, batch_file, assay="RNA", selection.method = "vst", nfeatures = 3000) {
  pca_dims = 1:npcs
  
  cat("NormalizeData ... \n")
  obj <- NormalizeData(obj, verbose = FALSE)
  
  cat("FindVariableFeatures ... \n")
  obj <- FindVariableFeatures(obj, selection.method = selection.method, nfeatures = nfeatures, verbose = FALSE)
  
  var.genes <- VariableFeatures(subobj)

  cat("ScaleData ... \n")
  all.genes <- rownames(obj)  
  obj <- ScaleData(obj, features = all.genes, verbose = FALSE)
  
  cat("RunPCA ... \n")
  obj <- RunPCA(object = obj, assay=assay, features = var.genes, verbose=FALSE)
  
  if(file.exists(batch_file)){
    cat("Setting batch ...\n")
    poolmap = get_batch_samples(batch_file, unique(obj$sample))
    obj$batch <- unlist(poolmap[obj$sample])
  }else if(!("batch" %in% colnames(obj@meta.data))){
    obj$batch <- obj$sample
  }
  
  cat("RunHarmony ... \n")
  obj <- RunHarmony(object = obj,
                    assay.use = assay,
                    reduction = "pca",
                    dims.use = pca_dims,
                    group.by.vars = "batch",
                    do_pca=FALSE)
  return(obj)
}

cluster_to_cell_type<-function(clusterDf){
  result=clusterDf[!duplicated(clusterDf$seurat_clusters),]
  result<-result[order(result$seurat_clusters),]
  result<-result[,colnames(result) != "cell"]
  rownames(result)<-result$seurat_clusters
  return(result)
}

sort_cell_type<-function(cts, sort_column){
  result<-cts[order(cts$seurat_clusters),]
  ct_uniq<-result[!duplicated(result[,sort_column]),]
  result[,sort_column]=factor(result[,sort_column], levels=ct_uniq[,sort_column])
  result<-result[order(result[,sort_column], result[,"seurat_clusters"]),]
  return(result)
}

plot_violin<-function(obj, features=c("FKBP1A", "CD79A")){
  library(reshape2)
  library(Seurat)
  
  glist=VlnPlot(all_obj, features=features, combine = F)
  
  gdata<-glist[[1]]$data
  for(idx in c(2:length(glist))){
    cdata=glist[[idx]]$data
    gdata<-cbind(gdata, cdata[,1])
    colnames(gdata)[ncol(gdata)]=colnames(cdata)[1]
  }
  
  mdata<-melt(gdata, id="ident")
  
  ggplot(mdata, aes(ident, value)) + 
    geom_violin(aes(fill = ident), trim=TRUE, scale="width") + 
    geom_jitter(width=0.5,size=0.5) + 
    facet_grid(variable~.) + 
    theme_classic() + NoLegend() + 
    xlab("Cluster") + ylab("Expression") +
    theme(strip.background=element_blank(),
          strip.text.y = element_text(angle =0))
}

preprocessing_rawobj<-function(rawobj, myoptions, prefix){
  Mtpattern= myoptions$Mtpattern
  rRNApattern=myoptions$rRNApattern
  Remove_rRNA<-ifelse(myoptions$Remove_rRNA == "0", FALSE, TRUE)
  Remove_MtRNA<-ifelse(myoptions$Remove_MtRNA == "0", FALSE, TRUE)
  
  nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
  nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
  nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
  mt_cutoff=as.numeric(myoptions$mt_cutoff)
  
  plot1 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_hline(yintercept = mt_cutoff, color="black")  + geom_vline(xintercept = nCount_cutoff, color="black")
  plot2 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  + geom_hline(yintercept = c( nFeature_cutoff_min, nFeature_cutoff_max), color="black")  + geom_vline(xintercept = nCount_cutoff, color="black") 
  p<-plot1+plot2
  
  png(paste0(prefix, ".qc.1.png"), width=3300, height=1500, res=300)
  print(p)
  dev.off()
  
  mt<-data.frame(mt=rawobj$percent.mt, Sample=rawobj$orig.ident, nFeature=log10(rawobj$nFeature_RNA), nCount=log10(rawobj$nCount_RNA))
  g1<-ggplot(mt, aes(y=mt,x=nCount) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_hline(yintercept = mt_cutoff, color="red")  + 
    geom_vline(xintercept = log10(nCount_cutoff), color="red") +
    ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
    facet_wrap(~Sample) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  
  png(paste0(prefix, ".qc.2.png"), width=3600, height=3000, res=300)
  print(g1)
  dev.off()
  
  g1<-ggplot(mt, aes(y=mt,x=nFeature) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_hline(yintercept = mt_cutoff, color="red")  + 
    geom_vline(xintercept = log10(nFeature_cutoff_min), color="red") +
    geom_vline(xintercept = log10(nFeature_cutoff_max), color="red") +
    ylab("Percentage of mitochondrial") + xlab("log10(number of feature)") +
    facet_wrap(~Sample) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  png(paste0(prefix, ".qc.3.png"), width=3600, height=3000, res=300)
  print(g1)
  dev.off()
  
  finalList<-list()
  
  #filter cells
  finalList$filter<-list(nFeature_cutoff_min=nFeature_cutoff_min,
                         nFeature_cutoff_max=nFeature_cutoff_max,
                         mt_cutoff=mt_cutoff,
                         nCount_cutoff=nCount_cutoff)
  
  rawCells<-data.frame(table(rawobj$orig.ident))
  
  rawobj<-subset(rawobj, subset = nFeature_RNA > nFeature_cutoff_min & nFeature_RNA<nFeature_cutoff_max & nCount_RNA > nCount_cutoff & percent.mt < mt_cutoff)
  
  if(Remove_rRNA){
    rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
    rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]
  }
  
  if(Remove_MtRNA){
    Mt.genes<- grep (pattern= Mtpattern,rownames(rawobj), value=TRUE )
    rawobj<-rawobj[!(rownames(rawobj) %in% Mt.genes),]
  }
  
  filteredCells<-data.frame(table(rawobj$orig.ident))
  qcsummary<-merge(rawCells, filteredCells, by = "Var1")
  colnames(qcsummary)<-c("Sample", "RawCell", "ValidCell")
  qcsummary$DiscardCell<-qcsummary$RawCell-qcsummary$ValidCell
  qcsummary$DiscardRate<-qcsummary$DiscardCell / qcsummary$RawCell
  write.csv(qcsummary, file=paste0(prefix, ".sample_cell.csv"), row.names=F)
  
  g<-VlnPlot(object = rawobj, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by="sample")
  png(paste0(prefix, ".qc.4.png"), width=3600, height=1600, res=300)
  print(g)
  dev.off()
  
  finalList$rawobj<-rawobj
  return(finalList)
}

output_integration_dimplot<-function(obj, outFile, has_batch_file){
  g<-FeaturePlot(obj, features="percent.mt") + ggtitle("Percentage of mitochondrial genes")
  width=2500
  ncol=1
  
  if("percent.hb" %in% colnames(obj@meta.data)){
    g2<-FeaturePlot(obj, features="percent.hb") + ggtitle("Percentage of hemoglobin genes")
    g<-g+g2
    width=width+2200
    ncol=ncol+1
  }

  if("percent.ribo" %in% colnames(obj@meta.data)){
    g3<-FeaturePlot(obj, features="percent.hb") + ggtitle("Percentage of ribosomal genes")
    g<-g+g3
    width=width+2200
    ncol=ncol+1
  }
  g=g+plot_layout(ncol=ncol)

  png(paste0(outFile, ".genes.png"), width=width, height=2000, res=300)
  print(g)
  dev.off()  
  
  mt<-data.frame(UMAP_1=obj@reductions$umap@cell.embeddings[,1], 
                UMAP_2=obj@reductions$umap@cell.embeddings[,2],
                Sample=obj$sample,
                Ident=obj$orig.ident)
  if(has_batch_file){
    if(!all(obj$sample == obj$batch)){
      mt$batch=obj$batch
    }else{
      has_batch_file=FALSE
    }
  }

  nSplit = length(unique(mt[,"Sample"]))
  nWidth=ceiling(sqrt(nSplit))
  nHeight=ceiling(nSplit / nWidth)
  width=nWidth * 600 + 200
  height=nHeight*600
  
  cat("draw pictures ... ")
  p<-draw_dimplot(mt, paste0(outFile, ".sample.png"), "Sample")
  if(!all(mt$Sample == mt$Ident)){
    p1<-draw_dimplot(mt, paste0(outFile, ".Ident.png"), "Ident")
    p<-p+p1
    width=width + nWidth * 600
  }

  if(has_batch_file){
    p2<-draw_dimplot(mt, paste0(outFile, ".batch.png"), "batch")
    p<-p+p2
    width=width+nWidth * 600
  }
  
  png(paste0(outFile, ".final.png"), width=width, height=height, res=300)
  print(p)
  dev.off()
}

read_bubble_genes<-function(bubble_file, allgenes=NA){
  library("readxl")
  
  genes <- read_xlsx(bubble_file, sheet = 1)
  for(idx in c(2:nrow(genes))){
    if(is.na(genes[idx,"Cell Type"])){
      genes[idx,"Cell Type"]=genes[idx-1,"Cell Type"]
    }
  }
  
  gene_names=genes$`Marker Gene`
  gene_names[gene_names=="PECAM"] = "PECAM1"
  gene_names[gene_names=="HGD1B"] = "HGD"
  gene_names[gene_names=="EpCAM"] = "EPCAM"
  gene_names[gene_names=="CD25"] = "IL2RA"
  gene_names[gene_names=="ACTAA2"] = "ACTA2"
  gene_names[gene_names=="MTND6"] = "MT-ND6"
  gene_names[gene_names=="FOXJ!"] = "FOXJ1"
  
  genes$`Marker Gene`<-gene_names
  
  miss_genes=genes$`Marker Gene`[!(genes$`Marker Gene` %in% allgenes)]
  writeLines(miss_genes, con="miss_gene.csv")
  
  if(!is.na(allgenes)){
    genes<-genes[genes$`Marker Gene` %in% allgenes,]
  }
  genes$`Cell Type`=factor(genes$`Cell Type`, levels=unique(genes$`Cell Type`))
  
  return(genes)
}

find_number_of_reduction<-function(obj, reduction="pca"){
  pct <- obj[[reduction]]@stdev / sum(obj[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  return(pcs)
}

get_seurat_average_expression<-function(SCLC, cluster_name){
  dd=GetAssayData(SCLC,slot="data")
  dobj=CreateSeuratObject(counts=dd)
  dobj$seurat_clusters=SCLC[[cluster_name]]
  result<-AverageExpression(dobj, slot="counts", group.by="seurat_clusters" )[[1]]
  rm(dd)
  rm(dobj)
  return(result)
}

get_bubble_plot<-function(obj, cur_res, cur_celltype, bubblemap_file, assay){
  allgenes=rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  gene_groups=split(genes_df$`Marker Gene`, genes_df$`Cell Type`)
  
  cell_type=obj@meta.data
  cell_type$cell_type <- cell_type[,cur_celltype]

  ct_levels<-c("B cells", "Plasma cells", "NK cells", "T cells", "Macrophages", "Dendritic cells", "Monocytes", "Mast cells", "Endothelial cells", "Fibroblasts", "Epithelial cells", "Basal cells", "Olfactory epithelial cells", "Ciliated cells")
  ct<-cell_type[!duplicated(cell_type$cell_type),]
  missed = ct$cell_type[!(ct$cell_type %in% ct_levels)]
  if(length(missed) > 0){
    ct_levels = c(ct_levels, missed)
  }
  ct_levels = ct_levels[ct_levels %in% ct$cell_type]
  cell_type$cell_type<-factor(cell_type$cell_type, levels=ct_levels)
  if(!is.na(cur_res)){
    cell_type<-cell_type[order(cell_type$cell_type, cell_type[,cur_res]),]
    cell_type$seurat_celltype_clusters=paste0(cell_type[,cur_res], " : ", cell_type$cell_type)
    cell_type$seurat_celltype_clusters=factor(cell_type$seurat_celltype_clusters, levels=unique(cell_type$seurat_celltype_clusters))
    group.by="seurat_celltype_clusters"
  }else{
    group.by="cell_type"
  }
  
  cell_type<-cell_type[colnames(obj),]
  obj@meta.data<-cell_type
  
  
  genes=unique(unlist(gene_groups))
  g<-DotPlot(obj, features=genes, assay=assay,group.by=group.by)
  gdata<-g$data

  data.plot<-NULL
  gn=names(gene_groups)[1]
  for(gn in names(gene_groups)){
    gs=gene_groups[[gn]]
    gdd<-gdata[gdata$features.plot %in% gs,]
    if(nrow(gdd)== 0){
      stop(gn)
    }
    gdd$feature.groups=gn
    data.plot<-rbind(data.plot, gdd)
  }
  
  data.plot$feature.groups=factor(data.plot$feature.groups, levels=names(gene_groups))
  
  color.by <- "avg.exp.scaled"
  scale.func <- scale_radius
  scale.min = NA
  scale.max = NA
  dot.scale = 6
  cols = c("lightgrey", "blue")
  
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = "Identity") +
    theme_cowplot() + 
    facet_grid(facets = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") + 
    theme(panel.spacing = unit(x = 1,units = "lines"), strip.background = element_blank()) + 
    scale_color_gradient(low = cols[1], high = cols[2])
  
  
  g=plot + 
    xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                                             axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                                             strip.background = element_blank(),
                                             strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
  return(g)
}

draw_bubble_plot<-function(obj, cur_res, cur_celltype, bubble_map_file, assay, prefix){
  g<-get_bubble_plot(obj, cur_res, cur_celltype, bubble_map_file, assay)
  png(paste0(prefix, ".bubblemap.png"), width=5500, height=3000, res=300)
  print(g)
  dev.off()
}
