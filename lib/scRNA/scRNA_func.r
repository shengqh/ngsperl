library(harmony)
library(cowplot)
library(Seurat)
library(tools)
library(scales)
library(ggplot2)
library(patchwork)
library(Matrix.utils)
library(parallel)

check_mc_cores<-function(mc.cores) {  
  if(.Platform$OS.type == "windows") {
    mc.cores=1
  }else{
    mc.cores=min(parallel::detectCores() - 1, max(1, mc.cores))
  }
  return(mc.cores)
}

to_numeric<-function(value, defaultValue){
  if(is.null(value)){
    return(defaultValue)
  }
  if(is.na(value)){
    return(defaultValue)
  }
  if(value == ""){
    return(defaultValue)
  }
  return(as.numeric(value))
}

is_one<-function(value, defaultValue=FALSE){
  if(is.null(value)){
    return(defaultValue)
  }
  if(is.na(value)){
    return(defaultValue)
  }
  return(value == '1')
}

is_file_empty<-function(filepath){
  if(is.null(filepath)){
    return(TRUE)
  }

  if(is.na(filepath)){
    return(TRUE)
  }

  if('' == filepath){
    return(TRUE)
  }

  if(!file.exists(filepath)){
    return(TRUE)
  }
  
  if(file.info(filepath)$size == 0){
    return(TRUE)
  }

  return(FALSE)
}

celltype_to_filename<-function(pct){
  return(gsub('[/:()?\ ]+', "_", pct))
}

get_hue_colors<-function(n, random_colors=TRUE, random.seed=20220606){
  ccolors<-hue_pal()(n)
  if(random_colors){
    x <- .Random.seed
    set.seed(random.seed)
    ccolors<-sample(ccolors, size=n)
    .Random.seed <- x
  }
  return(ccolors)
}

theme_rotate_x_axis_label <- function() {
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

theme_bw3 <- function (axis.x.rotate=F) { 
	result = theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      axis.line = element_line(colour = "black", size = 0.5)
    )
  if (axis.x.rotate){
    result = result + theme_rotate_x_axis_label()
  }

  return(result)
}

#https://github.com/satijalab/seurat/issues/1836
#For visualization, using sctransform data is also fine.

MyDoHeatMap<-function(object, max_cell=5000, ...){
  if(ncol(obj) > max_cell){
    subsampled <- obj[, sample(colnames(obj), size=max_cell, replace=F)]
    g<-DoHeatmap(subsampled, ...)
  }else{
    g<-DoHeatmap(object, ...)
  }
  return(g)
}

MyFeaturePlot<-function(object, assay="RNA", ...){
  old_assay=DefaultAssay(object)
  DefaultAssay(object)=assay
  g=FeaturePlot(object=object, ...)
  DefaultAssay(object)=old_assay
  g
}

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

get_assay<-function(by_sctransform, by_integration, by_harmony){
  assay=ifelse(by_sctransform, "SCT", "RNA")
  if(by_integration & (!by_harmony)){
    assay="integrated"
  }
  return(assay)
}

calc_weight<-function(cellType){
  freq<-sort((table(unlist(cellType)))/length(cellType))
  weight<-1+sqrt((max(freq)-freq)/(max(freq)-min(freq)))
  return(weight)
}

get_celltype_map<-function( root_celltypes, HLA_panglao5_file ){
  cts<-read.table(HLA_panglao5_file, sep="\t", header=T, row.names=1)
  cmap=cts$Layer4
  names(cmap) = rownames(cts)

  layer="Layer3"
  for(layer in c("Layer1", "Layer2", "Layer3")){
    lcts<-cts[cts[,layer] %in% root_celltypes,]
    cmap[rownames(lcts)]<-lcts[,layer]
  }

  return(cmap)
}

# 20221211, remove remove_subtype_str and HLA_panglao5_file parameters
read_cell_markers_file<-function(panglao5_file, species, curated_markers_file=NULL){
  #preparing cell activity database
  marker<-data.frame(fread(panglao5_file))

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
  
  if(!is_file_empty(curated_markers_file)){
    curated_markers_df<-read.table(curated_markers_file, sep="\t", header=F, stringsAsFactors=F)
    curated_markers_celltype<-split(curated_markers_df$V2, curated_markers_df$V1)
    for(cmct in names(curated_markers_celltype)){
      cellType[[cmct]]=curated_markers_celltype[[cmct]]
    }
    weight=calc_weight(cellType)
  } 

  return(list(cellType=cellType, weight=weight))
}

update_tiers<-function(tiers, celltype_map){
  rownames(tiers)<-tiers$Celltype.name

  update_layers=paste0("Layer", c(1,2,3,4))
  tiers[names(celltype_map),update_layers]<-tiers[celltype_map,update_layers]

  return(tiers)
}

read_tiers<-function(HLA_panglao5_file, remove_subtype_str = NULL){
  tiers<-read.table(HLA_panglao5_file, sep="\t", header=T)
  if(!is.null(remove_subtype_str)) {
    remove_subtype_list<-unique(unlist(strsplit(remove_subtype_str, ",")))
    if(length(remove_subtype_list) > 0){
      celltype_map = get_celltype_map(remove_subtype_list, HLA_panglao5_file)
      celltype_map = celltype_map[names(celltype_map) != celltype_map]
      tiers = update_tiers(tiers, celltype_map)
    }
  }

  return(tiers)
}

get_summary_layer<-function(tiers, layer, remove_subtype_str=NULL){
  result<-unlist(split(tiers[,layer], tiers$Celltype.name))

  if(!is.null(remove_subtype_str)){
    remove_subtype_list<-unique(unlist(strsplit(remove_subtype_str, ",")))
    result[remove_subtype_list]<-remove_subtype_list
  }

  return(result)
}

init_celltype_markers<-function(panglao5_file, species, curated_markers_file=NULL, HLA_panglao5_file, layer="Layer4", remove_subtype_str=NULL, combined_celltype_file=NULL){
  tiers<-read_tiers(HLA_panglao5_file, remove_subtype_str)

  cell_activity_database<-read_cell_markers_file(panglao5_file=panglao5_file, 
                                                 species=species, 
                                                 curated_markers_file=curated_markers_file)

  cell_activity_database$cellType=cell_activity_database$cellType[names(cell_activity_database$cellType) %in% tiers$Celltype.name]
  tiers<-tiers[tiers$Celltype.name %in% names(cell_activity_database$cellType),]

  celltype_map=get_summary_layer(tiers, layer, remove_subtype_str)

  combined_ct=NA
  combined_ct_source<-NA
  if(!is_file_empty(combined_celltype_file)){
    cts<-read.table(combined_celltype_file, header=F, sep="\t", stringsAsFactors = F)
    combined_ct<-unlist(split(cts$V2, cts$V1))

    layer2map<-celltype_map
    combined_ct_in_layer<-unique(combined_ct[combined_ct %in% names(layer2map)])
    combined_ct[combined_ct_in_layer]<-combined_ct_in_layer
    layer2map[layer2map %in% names(combined_ct)]<-combined_ct[layer2map[layer2map %in% names(combined_ct)]]
    combined_ct_source<-split(names(combined_ct), combined_ct)
  }

  return(list(tiers=tiers, cell_activity_database=cell_activity_database, celltype_map=celltype_map, combined_celltypes=combined_ct, combined_celltype_source=combined_ct_source))
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
  markers=FindMarkers(object, assay="RNA", ident.1=ident.1, ident.2=ident.2, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  markers=markers[markers$p_val_adj < 0.01,]
  return(markers)
}

find_all_markers<-function(object, by_sctransform, min.pct = 0.5, logfc.threshold = 0.6){
  markers=FindAllMarkers(object, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
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

run_cluster<-function(object, pca_dims, resolution, random.seed, reduction="pca", min.pct = 0.5, logfc.threshold = 0.6){
  object=run_cluster_only(object, pca_dims, resolution, random.seed, reduction)
  markers <- FindAllMarkers(object, assay="RNA", only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)
  markers <- markers[markers$p_val_adj < 0.05,]
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

do_normalization<-function(obj, selection.method, nfeatures, vars.to.regress, scale.all=FALSE, essential_genes=NULL) {
  DefaultAssay(obj)<-"RNA"

  cat("NormalizeData ... \n")
  obj <- NormalizeData(obj, verbose = FALSE)
  cat("NormalizeData done ... \n")
  
  cat("FindVariableFeatures ... \n")
  obj <- FindVariableFeatures(obj, selection.method = selection.method, nfeatures = nfeatures, verbose = FALSE)
  cat("FindVariableFeatures done ... \n")
  
  if(scale.all){
    features = rownames(obj@assays$RNA@data)
  }else{
    features=VariableFeatures(obj)
    if (!is.null(essential_genes)){
      features = unique(c(features, essential_genes))
    }
  }
  cat("ScaleData: ", length(features), "genes ...\n")
  obj <- ScaleData(obj, vars.to.regress=vars.to.regress, features=features, verbose = FALSE)
  cat("ScaleData done ... \n")

  stopifnot(dim(obj@assays$RNA@scale.data)[1] > 0)

  return(obj)
}

do_sctransform<-function(rawobj, vars.to.regress, return.only.var.genes=FALSE, mc.cores=1) {
  mc.cores = check_mc_cores(mc.cores)

  print("performing SCTransform ...")
  nsamples=length(unique(rawobj$orig.ident))
  if(nsamples > 1){
    print("  split objects ...")
    objs<-SplitObject(object = rawobj, split.by = "orig.ident")

    print("  perform sctransform ...")
    objs<-mclapply(objs, function(x){
      print(paste0("    sctransform ", unique(x$orig.ident), " ..."))
      x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE)
      return(x)
    },mc.cores=mc.cores)  
    print("  sctransform done")

    print("  merge samples ...")
    obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
    #https://github.com/satijalab/seurat/issues/2814
    VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data)
    rm(objs)
    return(obj)
  }else{
    print("  perform sctransform ...")
    rawobj<-SCTransform(rawobj, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE)
    print("  sctransform done")
    return(rawobj)
  }
}

do_harmony<-function(obj, by_sctransform, vars.to.regress, has_batch_file, batch_file, pca_dims, essential_genes=NULL, mc.cores=1){
  if(by_sctransform){
    #now perform sctranform
    obj<-do_sctransform(obj, vars.to.regress=vars.to.regress,mc.cores=mc.cores)
    assay="SCT"
  }else{
    assay="RNA"
  }

  #no matter if we will use sctransform, we need normalized RNA assay for visualization and cell type annotation
  #data slot for featureplot, dotplot, cell type annotation and scale.data slot for heatmap
  #we need to do sctransform first, then do RNA assay normalization and scale, otherwise, after sctransform, the scale.data slot will be discarded.
  obj<-do_normalization(obj, selection.method="vst", nfeatures=3000, vars.to.regress=vars.to.regress, scale.all=FALSE, essential_genes=essential_genes)

  DefaultAssay(obj)<-assay

  print("RunPCA ...")
  obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)

  if(has_batch_file){
    print("Setting batch ...")
    poolmap = get_batch_samples(batch_file, unique(obj$sample))
    obj$batch <- unlist(poolmap[obj$sample])
  }else if(!("batch" %in% colnames(obj@meta.data))){
    obj$batch <- obj$sample
  }

  print("RunHarmony ...")
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
  
  rawCells<-data.frame(table(rawobj$orig.ident))
  
  if(Remove_rRNA){
    rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
    rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]
  }

  if(Remove_MtRNA){
    Mt.genes<- grep (pattern= Mtpattern,rownames(rawobj), value=TRUE )
    rawobj<-rawobj[!(rownames(rawobj) %in% Mt.genes),]
  }
  
  plot1 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
    geom_hline(yintercept = mt_cutoff, color="black")  + 
    geom_vline(xintercept = nCount_cutoff, color="black") +
    scale_y_continuous(breaks = seq(0, 100, by = 10))
  plot2 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")  + 
    geom_hline(yintercept = c( nFeature_cutoff_min, nFeature_cutoff_max), color="black")  + 
    geom_vline(xintercept = nCount_cutoff, color="black") 
  p<-plot1+plot2
  
  png(paste0(prefix, ".qc.1.png"), width=3300, height=1500, res=300)
  print(p)
  dev.off()

  mt<-data.frame(mt=rawobj$percent.mt, Sample=rawobj$orig.ident, nFeature=log10(rawobj$nFeature_RNA), nCount=log10(rawobj$nCount_RNA))

  nsample<-length(unique(mt$Sample))
  ncol=ceiling(sqrt(nsample))
  nrow=ceiling(nsample/ncol)
  width=min(10000, ncol * 1200)
  height=min(10000, nrow * 1000)

  g1<-ggplot(mt, aes(y=mt,x=nCount) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_hline(yintercept = mt_cutoff, color="red")  + 
    geom_vline(xintercept = log10(nCount_cutoff), color="red") +
    ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    facet_wrap(~Sample) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  png(paste0(prefix, ".qc.2.png"), width=width, height=height, res=300)
  print(g1)
  dev.off()
  
  g1<-ggplot(mt, aes(y=mt,x=nFeature) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_hline(yintercept = mt_cutoff, color="red")  + 
    geom_vline(xintercept = log10(nFeature_cutoff_min), color="red") +
    geom_vline(xintercept = log10(nFeature_cutoff_max), color="red") +
    ylab("Percentage of mitochondrial") + xlab("log10(number of feature)") +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    facet_wrap(~Sample) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  png(paste0(prefix, ".qc.3.png"), width=width, height=height, res=300)
  print(g1)
  dev.off()

  finalList<-list()
  
  #filter cells
  finalList$filter<-list(nFeature_cutoff_min=nFeature_cutoff_min,
                         nFeature_cutoff_max=nFeature_cutoff_max,
                         mt_cutoff=mt_cutoff,
                         nCount_cutoff=nCount_cutoff,
                         Remove_rRNA=Remove_rRNA,
                         Remove_MtRNA=Remove_MtRNA)

  rawobj<-subset(rawobj, subset = nFeature_RNA >= nFeature_cutoff_min & nFeature_RNA <= nFeature_cutoff_max & nCount_RNA >= nCount_cutoff & percent.mt <= mt_cutoff)

  filteredCells<-data.frame(table(rawobj$orig.ident))
  qcsummary<-merge(rawCells, filteredCells, by = "Var1")
  colnames(qcsummary)<-c("Sample", "RawCell", "ValidCell")
  qcsummary$DiscardCell<-qcsummary$RawCell-qcsummary$ValidCell
  qcsummary$DiscardRate<-qcsummary$DiscardCell / qcsummary$RawCell
  write.csv(qcsummary, file=paste0(prefix, ".filtered.cell.csv"), row.names=F)
  
  g<-VlnPlot(object = rawobj, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by="orig.ident")
  png(paste0(prefix, ".qc.4.png"), width=3600, height=1600, res=300)
  print(g)
  dev.off()
  
  finalList$rawobj<-rawobj
  return(finalList)
}

output_integration_dimplot<-function(obj, outFile, has_batch_file, qc_genes=NULL){
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
    g3<-FeaturePlot(obj, features="percent.ribo") + ggtitle("Percentage of ribosomal genes")
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
    if(!all(obj$orig.ident == obj$batch)){
      mt$batch=obj$batch
    }else{
      has_batch_file=FALSE
    }
  }

  nSplit = length(unique(mt[,"Ident"]))
  nWidth=ceiling(sqrt(nSplit))
  nHeight=ceiling(nSplit / nWidth)
  width=nWidth * 600 + 200
  height=nHeight*600
  
  cat("draw pictures ... ")
  draw_feature_qc(paste0(outFile, ".Ident"), obj, "orig.ident")

  p<-draw_dimplot(mt, paste0(outFile, ".Ident.png"), "Ident")
  if(!all(mt$Sample == mt$Ident)){
    draw_feature_qc(paste0(outFile, ".sample"), obj, "sample")
    p1<-draw_dimplot(mt, paste0(outFile, ".sample.png"), "Sample")
    p<-p+p1
    width=width + nWidth * 600
  }

  if(has_batch_file){
    draw_feature_qc(paste0(outFile, ".batch"), obj, "batch")
    p2<-draw_dimplot(mt, paste0(outFile, ".batch.png"), "batch")
    p<-p+p2
    width=width+nWidth * 600
  }
  
  png(paste0(outFile, ".final.png"), width=width, height=height, res=300)
  print(p)
  dev.off()

  if(!is.null(qc_genes) & qc_genes != ''){
    genes<-unlist(strsplit( qc_genes, ',' ))
    g<-FeaturePlot(obj, genes, split.by="orig.ident")
    png(paste0(outFile, ".qc_genes.png"), width=3000, height=6000, res=300)
    print(g)
    dev.off()
  }
}

read_bubble_genes<-function(bubblemap_file, allgenes=c()){
  library("readxl")
  library("tidyr")
  
  genes <- read_xlsx(bubblemap_file, sheet = 1)
  colnames(genes)[colnames(genes) == "Marker Gene"] = "gene"
  colnames(genes)[colnames(genes) == "Cell Type"] = "cell_type"

  for(idx in c(2:nrow(genes))){
    if(is.na(genes[idx,"cell_type"])){
      genes[idx,"cell_type"]=genes[idx-1,"cell_type"]
    }
  }
  
  genes<-separate_rows(genes, gene)
  
  gene_names=genes$gene
  gene_names[gene_names=="PECAM"] = "PECAM1"
  gene_names[gene_names=="HGD1B"] = "HGD"
  gene_names[gene_names=="EpCAM"] = "EPCAM"
  gene_names[gene_names=="CD25"] = "IL2RA"
  gene_names[gene_names=="ACTAA2"] = "ACTA2"
  gene_names[gene_names=="MTND6"] = "MT-ND6"
  gene_names[gene_names=="FOXJ!"] = "FOXJ1"
  
  genes$gene<-gene_names
  
  if(length(allgenes) > 0){
    miss_genes=setdiff(genes$gene, allgenes)
    writeLines(miss_genes, con="miss_gene.csv")

    genes<-genes[genes$gene %in% allgenes,]
  }
  genes$cell_type=factor(genes$cell_type, levels=unique(genes$cell_type))
  
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

get_seurat_average_expression<-function(SCLC, cluster_name, assay="RNA"){
  dd=GetAssayData(SCLC, assay=assay, slot="data")
  dobj=CreateSeuratObject(counts=dd)
  dobj$seurat_clusters=SCLC[[cluster_name]]
  result<-AverageExpression(dobj, slot="counts", group.by="seurat_clusters" )[[1]]

  cnames = unique(dobj$seurat_clusters)
  cnames = cnames[order(cnames)]
  colnames(result) <- cnames

  rm(dd)
  rm(dobj)
  return(result)
}

get_dot_plot<-function(obj, group.by, gene_groups, assay="RNA", rotate.title=TRUE ){
  genes=unique(unlist(gene_groups))
  g<-DotPlot(obj, features=genes, assay=assay, group.by=group.by)
  
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
  
  
  g=plot + xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                             axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                                             strip.background = element_blank())
  
  if(rotate.title){
    g=g+theme(strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
  }
  return(g)
}

get_bubble_plot<-function(obj, cur_res, cur_celltype, bubblemap_file, assay="RNA", orderby_cluster=FALSE, split.by=NULL, rotate.title=TRUE, group.by=NULL, use_blue_yellow_red=TRUE){
  allgenes=rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  gene_groups=split(genes_df$gene, genes_df$cell_type)

  if(is.null(group.by)){
    if(is.null(cur_celltype)){
      stop("cur_celltype cannot be null in get_bubble_plot if group.by is null.")
    }
    if(is.null(cur_res) || is.na(cur_res)){
      cur_res = paste0(cur_celltype, "_cluster")
      obj<-build_dummy_cluster(obj, label.by=cur_celltype, new_cluster_name=cur_res)
    }

    cell_type=obj@meta.data
    cell_type$cell_type <- cell_type[,cur_celltype]

    # ct_levels<-c("B cells", "Plasma cells", "NK cells", "T cells", "Macrophages", "Dendritic cells", "Monocytes", "Mast cells", "Endothelial cells", "Fibroblasts", "Epithelial cells", "Basal cells", "Olfactory epithelial cells", "Ciliated cells")
    # ct<-cell_type[!duplicated(cell_type$cell_type),]
    # missed = ct$cell_type[!(ct$cell_type %in% ct_levels)]
    # if(length(missed) > 0){
    #   ct_levels = c(ct_levels, as.character(missed))
    # }
    # ct_levels = ct_levels[ct_levels %in% ct$cell_type]
    # cell_type$cell_type<-factor(cell_type$cell_type, levels=ct_levels)

    if (! is.null(split.by)){ 
      if(orderby_cluster){
        cell_type<-cell_type[order(cell_type[,cur_res], cell_type[,split.by]),]
      }else{
        cell_type<-cell_type[order(cell_type$cell_type, cell_type[,cur_res], cell_type[,split.by]),]
      }
      cell_type$seurat_celltype_clusters=paste0(cell_type[,cur_res], ": ", cell_type$cell_type, ": ", cell_type[,split.by])
    }else{
      if(orderby_cluster){
        cell_type<-cell_type[order(cell_type[,cur_res]),]
      }else{
        cell_type<-cell_type[order(cell_type$cell_type, cell_type[,cur_res]),]
      }
      cell_type$seurat_celltype_clusters=paste0(cell_type[,cur_res], ": ", cell_type$cell_type)
    }
    cell_type$seurat_celltype_clusters=factor(cell_type$seurat_celltype_clusters, levels=unique(cell_type$seurat_celltype_clusters))
    group.by="seurat_celltype_clusters"
    
    cell_type<-cell_type[colnames(obj),]
    obj@meta.data<-cell_type
  }

  g<-get_dot_plot(obj, group.by, gene_groups, assay, rotate.title=rotate.title)
  
  if(use_blue_yellow_red){
    g <- g + scale_color_gradient2(low="blue", mid="yellow", high="red")
  }
  
  return(g)
}

draw_bubble_plot<-function(obj, cur_res, cur_celltype, bubble_map_file, prefix, width=5500, height=3000, rotate.title=TRUE){
  g<-get_bubble_plot(obj, cur_res, cur_celltype, bubble_map_file, rotate.title=rotate.title)
  png(paste0(prefix, ".bubblemap.png"), width=width, height=height, res=300)
  print(g)
  dev.off()
}

get_celltype_markers<-function(medianexp,cellType,weight,combined_ct_source=NA,layer_map=NA){
  exp_z<-scale(medianexp)
  genenames<-rownames(medianexp)   
  ctnames<-colnames(medianexp)
  j=7
  
  res<-list()
  for (j in 1: dim(medianexp)[2]){
    clusterexp<-medianexp[,j] 
    clusterexp_z<-exp_z[,j]
    
    if(all(is.na(combined_ct_source))){
      orig_ctnames=ctnames[j]
    }else{
      orig_ctnames=unique(c(ctnames[j], combined_ct_source[[ctnames[j]]]))
    }

    valid_ctnames=orig_ctnames[orig_ctnames %in% names(cellType)]
    if(length(valid_ctnames) == 0){
      if(!all(is.na(layer_map))){
        valid_ctnames=names(layer_map)[layer_map%in%orig_ctnames]
      }
    }

    if(length(valid_ctnames) == 0){
      next
    }
    
    ctgenes<-unique(unlist(cellType[valid_ctnames]))
    stopifnot(length(ctgenes) > 0)
    
    weight_ss<-weight[names(weight)%in%ctgenes]
    ind<-match(names(weight_ss),genenames)
    exp_ss<-clusterexp_z[ind[!is.na(ind)]]
    weight_ss<-weight_ss[!is.na(ind)]
    
    weighted_exp<-exp_ss*weight_ss
    weighted_exp<-weighted_exp[order(weighted_exp, decreasing = T)]
    res[ctnames[j]] = list(weighted_exp)
  }
  
  return(res)
}

get_celltype_marker_bubble_plot<-function(obj, group.by, cellType, weight, n_markers=5, combined_ct_source=NA, layer_map=NA) {
  medianexp=get_seurat_average_expression(obj, group.by)
  medianexp<-medianexp[rowSums(medianexp)>0,]
  #medianexp<-medianexp[,colnames(medianexp)[order(colnames(medianexp))]]
  
  markers<-get_celltype_markers(medianexp, cellType, weight, combined_ct_source=combined_ct_source, layer_map=layer_map)
  
  while(TRUE){
    gene_groups<-lapply(markers, function(x){
      names(x[1:n_markers])
    })

    gc=table(unlist(gene_groups))
    gc=gc[gc>1]
    #print(gc)
    if(length(gc) > 0){
      markers<-lapply(markers, function(x){
        x[!(names(x) %in% names(gc))]
      })
    }else{
      break
    }
  }

  g<-get_dot_plot(obj, group.by, gene_groups)
  return(g)
}

find_best_resolution<-function(subobj, clusters, min.pct, logfc.threshold, min_markers) {  
  cluster = clusters[1]
  lastCluster = clusters[1]
  for(cluster in clusters){
    cat("  ", cluster, "\n")
    Idents(subobj)<-cluster
    if(length(unique(Idents(subobj))) == 1){
      lastCluster = cluster
      next
    }
    
    markers=FindAllMarkers(subobj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    markers=markers[markers$p_val_adj < 0.05,]
    nmarkers=unlist(lapply(unique(Idents(subobj)), function(x){sum(markers$cluster==x)}))
    if(all(nmarkers >= min_markers)){
      lastCluster = cluster
    }else{
      break
    }
  }
  return(lastCluster)
}
  
read_object<-function(obj_file, meta_rds=NULL, columns=NULL){
  obj=readRDS(obj_file)
  if(is.list(obj)){
    obj<-obj$obj
  }
  
  if(!is.null(meta_rds)){
    if (meta_rds != ""){
      if(file_ext(meta_rds) == "rds"){
        meta.data=readRDS(meta_rds)
        if(any(colnames(obj) != rownames(meta.data))){
          obj=subset(obj, cells=rownames(meta.data))
        }
        if(all(colnames(obj@meta.data) %in% colnames(meta.data))){
          obj@meta.data = meta.data
          return(obj)
        }
      }else{
        meta.data=read.csv(meta_rds, stringsAsFactors = F, row.names=1)
      }
      if(is.null(columns)){
        columns = colnames(meta.data)
      }
      for(column in columns){
        obj=AddMetaData(obj, meta.data[,column], col.name=column)
      }
    }
  }
  return(obj)
}

output_ElbowPlot<-function(obj, outFile, reduction){
  p<-ElbowPlot(obj, ndims = 40, reduction = reduction)
  png(paste0(outFile, ".elbowplot.", reduction, ".png"), width=1500, height=1200, res=300)
  print(p)
  dev.off()
}

draw_feature_qc<-function(prefix, rawobj, ident_name) {
  Idents(rawobj)<-ident_name

  nsample<-length(unique(unlist(rawobj[[ident_name]])))
  
  feats<-c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo")
  if("percent.hb" %in% colnames(rawobj@meta.data)){
    feats<-c(feats, "percent.hb")
  }

  png(file=paste0(prefix, ".qc.violin.png"), width=6000, height=4000, res=300)
  g<-VlnPlot(rawobj, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
  print(g)
  dev.off()

  if('umap' %in% names(rawobj@reductions)){
    nfeature<-length(feats)

    by.col=nfeature>=nsample
    g<-FeaturePlot(rawobj, feats, split.by=ident_name, reduction="umap", order=T, by.col=by.col)
    if(by.col){
      width = nsample * 800
      height = nfeature * 700
    }else{
      width = nfeature * 800
      height = nsample * 700
    }

    png(paste0(prefix, ".qc.exp.png"), width=width, height=height, res=300)
    print(g)
    dev.off()  
  }

  png(file=paste0(prefix, ".qc.png"), width=2600, height=1200, res=300)
  p1 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
  p2 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
  p<-p1+p2+plot_layout(ncol=2)
  print(p)
  dev.off()
  
  mt<-data.frame(mt=rawobj$percent.mt, Sample=unlist(rawobj[[ident_name]]), nFeature=log10(rawobj$nFeature_RNA), nCount=log10(rawobj$nCount_RNA))
  nwidth=ceiling(sqrt(nsample))
  nheight=ceiling(nsample/nwidth)
  
  png(file=paste0(prefix, ".qc.read.png"), width=min(20000, 500 * nwidth + 300), height=min(10000, 500*nheight), res=300)
  p1<-ggplot(mt, aes(y=mt,x=nCount) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
    facet_wrap(Sample~.) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  print(p1)
  dev.off()

  png(file=paste0(prefix, ".qc.feature.png"), width=min(20000, 500 * nwidth + 300), height=min(10000, 500*nheight), res=300)
  p2<-ggplot(mt, aes(y=mt,x=nFeature) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ylab("Percentage of mitochondrial") + xlab("log10(number of feature)") +
    facet_wrap(Sample~.) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  print(p2)
  dev.off()

  ct<-as.data.frame(table(rawobj[[ident_name]]))
  colnames(ct)<-c("Sample","Cell")

  if("project" %in% colnames(rawobj@meta.data)){
    is_hto=any(rawobj$project != rawobj$orig.ident)
  }else{
    is_hto=FALSE
  }

  if("sample" %in% colnames(rawobj@meta.data)){
    is_merged=any(rawobj$sample != rawobj$orig.ident)
  }else{
    is_merged=FALSE
  }

  if(ident_name == "orig.ident" & is_hto){
    meta<-rawobj@meta.data
    meta<-meta[!duplicated(meta$orig.ident),,drop=F]
    smap=split(as.character(meta$project), as.character(meta$orig.ident))
    ct$Set=unlist(smap[ct$Sample])
  }

  if(ident_name == "sample" & is_merged){
    meta<-rawobj@meta.data
    meta<-meta[!duplicated(meta$sample),,drop=F]
    smap=split(as.character(meta$orig.ident), as.character(meta$sample))
    ct$Set=unlist(smap[ct$Sample])
  }
  write.table(ct, paste0(prefix, ".cell.txt"), sep="\t", row.names=F)
  
  if("Set" %in% colnames(ct)){
    g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Set))
  }else{
    g<-ggplot(ct, aes(x=Sample, y=Cell))
  }
  g<-g + geom_bar(stat="identity") + theme_bw3(axis.x.rotate = T)
  png(paste0(prefix, ".cell.bar.png"), width=max(3000, nrow(ct) * 100), height=2000, res=300)
  print(g)
  dev.off()
}

myScaleData<-function(object, features, assay, ...){
  scaled.genes<-rownames(object[[assay]]@scale.data)
  if(!all(features %in% scaled.genes)){
    new.genes<-unique(features, scaled.genes)
    object=ScaleData(object, features=new.genes, assay=assay, ... )
  }
  return(object)
}

get_top10_markers<-function(markers){
  markers=markers[markers$p_val_adj < 0.05,]
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
  return(top10)
}

factorize_layer<-function(obj, layer){
  ldata<-unlist(obj[[layer]])
  lt=table(ldata)
  lt<-lt[order(lt, decreasing=T)]
  levels=names(lt)
  obj[[layer]]=factor(ldata, levels=levels)
  return(obj)
}
    
unfactorize_layer<-function(obj, layer){
  obj[[layer]]=as.character(unlist(obj[[layer]]))
  return(obj)
}

get_groups_dot<-function(subobj, group1, group2){
  ct<-table(unlist(subobj[[group1]]))
  ct_names<-paste0(names(ct), " (", ct, ")")
  names(ct_names)<-names(ct)
  tbl<-table(unlist(subobj[[group1]]), unlist(subobj[[group2]]))
  mtbl<-reshape2::melt(tbl)
  colnames(mtbl)<-c(group1, group2, "Cell")
  mtbl[,group1]<-as.character(mtbl[,group1])
  mtbl[,group1]<-factor(ct_names[mtbl[,group1]], levels=ct_names)
  g<-ggplot(mtbl, aes_string(group1, "Cell", fill=group2)) + geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + xlab("") + ylab("")
  return(g)
}

RunMultipleUMAP<-function(subobj, nn=c(30,20,10), min.dist=c(0.3,0.1,0.05), curreduction="PCA", cur_pca_dims=c(1:30)){
  ops=data.frame("nn"=nn, "min.dist"=min.dist)
  umap_names<-c()
  for(idx in c(1:nrow(ops))){
    nn=ops[idx, 'nn']
    #nn<-min(nn, u_n_neighbors)
    min.dist=ops[idx, 'min.dist']
    umap_name = paste0("umap_nn", nn, "_dist", min.dist)
    cat(umap_name, "\n")
    umap_names<-c(umap_names, umap_name)
    umap_key = paste0("UMAPnn", nn, "dist", min.dist * 100, "_")
    subobj<-RunUMAP(object = subobj, reduction=curreduction, reduction.key=umap_key, reduction.name=umap_name, n.neighbors=nn, min.dist=min.dist, dims=cur_pca_dims, verbose = FALSE)
  }
  return(list(obj=subobj, umap_names=umap_names))
}

get_dim_plot<-function(obj, group.by, label.by, label=T, title=label.by, legend.title=label.by, reduction="umap", split.by=NULL, ncol=1, random_colors=TRUE, ...){
  labels<-obj@meta.data[,c(group.by, label.by)]
  labels<-labels[!duplicated(labels[,group.by]),]
  labels<-labels[order(labels[,group.by]),]
  cts<-as.character(labels[,label.by])

  ngroups<-length(unlist(unique(obj[[group.by]])))
  scolors = get_hue_colors(ngroups, random_colors)

  g<-DimPlot(obj, group.by=group.by, label=label, reduction=reduction, split.by=split.by, ...)+ 
    scale_color_manual(legend.title, values=scolors, labels = cts, guide = guide_legend(ncol=ncol)) + ggtitle(title)
  return(g)
}

build_dummy_cluster<-function(obj, label.by, new_cluster_name, new_cluster_name_label=paste0(new_cluster_name, "_label")){
  groups<-as.character(obj@meta.data[,label.by])
  gt<-table(groups)
  gt<-gt[order(gt, decreasing=T)]
  dummy_cluster<-c(0:(length(gt)-1))
  names(dummy_cluster)<-names(gt)
  dc<-factor(dummy_cluster[groups], levels=dummy_cluster)
  obj[[new_cluster_name]]<-dc
  obj[[new_cluster_name_label]]<-paste0(obj@meta.data[,new_cluster_name], ": ", groups)
  return(obj)
}

get_dim_plot_labelby<-function(obj, label.by, title=label.by, label=T, legend.title=label.by, reduction="umap", split.by=NULL, ncol=1, label_has_cluster=FALSE){
  group.by="dummy_cluster"
  group.label="dummy_label"

  obj<-build_dummy_cluster(obj, label.by, group.by, group.label)

  group.label=ifelse(label_has_cluster, label.by, group.label)

  g<-get_dim_plot(obj, group.by=group.by, label.by=group.label, label=label, title=title, legend.title=legend.title, reduction=reduction, split.by=split.by, ncol=ncol)
  return(g)
}

get_highlight_cell_plot<-function(obj, group.by, reduction="umap", reorder=TRUE) {
  cts<-table(obj[[group.by]])
  if(reorder){
    cts<-cts[order(cts, decreasing = T)]
  }
  
  g<-NULL
  for (ct in names(cts)){
    ct_count<-cts[ct]
    pct<-paste0(ct, "(", ct_count, ")")
    cells<-colnames(obj)[obj[[group.by]] == ct]
    g0<-DimPlot(obj, label=F, cells.highlight =cells, reduction = reduction) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", ct)) + NoLegend()
    if(is.null(g)){
      g<-g0
    }else{
      g<-g+g0
    }
  }
  
  return(list(g=g, cts=cts))
}

add_x_y_axis<-function(g){
  g<-g+annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  return(g)
}

save_highlight_cell_plot<-function(filename, obj, group.by, reduction="umap", reorder=TRUE, title=""){
  res<-get_highlight_cell_plot(obj, group.by = group.by, reduction = reduction, reorder=reorder)
  g<-res$g

  if(title != ""){
    g<-g+plot_annotation(title=title, theme = theme(plot.title = element_text(hjust = 0.5, size = 16)))
  }
  cts<-res$cts
  
  ncol<-ceiling(sqrt(length(cts)))
  nrow<-ceiling(length(cts)/ncol)
  
  width=1500 * ncol
  height=1500 * nrow
  
  g<-g+plot_layout(ncol=ncol)
  png(filename, width=width, height=height, res=300)
  print(g)
  dev.off()
}

sumcount<-function(ct_count, groupings){
  rescount<-t(aggregate.Matrix(t(ct_count), groupings=groupings, fun="sum"))
  psum<-rowSums(rescount)
  rescount<-rescount[psum>0,]
  return(rescount)
}

get_seurat_sum_count<-function(obj, cluster_name, min_cell_per_sample=1){
  clusterDf<-obj@meta.data
  cts = as.character(unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name]))
  prefixList<-celltype_to_filename(cts)

  res_files=c()

  idx<-1
  for (idx in c(1:length(cts))){
    ct = cts[idx]
    cat(ct, "\n")
    
    prefix = prefixList[idx]
    
    clusterCt<-clusterDf[clusterDf[,cluster_name] == ct,]
    de_obj<-subset(obj, cells=rownames(clusterCt))

    dtb<-table(de_obj$orig.ident)

    snames<-names(dtb)[dtb < min_cell_per_sample]
    if(length(snames) > 0){
      cat("those samples were excluded due to cell less than", min_cell_per_sample, ": ", paste(snames, collapse = ","), "\n")
      dtb<-dtb[dtb >= min_cell_per_sample]
      de_obj<-subset(de_obj, orig.ident %in% names(dtb))
    }

    ct_count<-de_obj@assays$RNA@counts
    groupings<-unlist(de_obj$orig.ident)
    p_count<-sumcount(ct_count, groupings)

    p_file=paste0(prefix, ".pseudo_count.csv")
    
    write.csv(p_count, p_file)
    
    res_files<-c(res_files, file_path_as_absolute(p_file))
  }

  res_df<-data.frame("cluster"=cts, "prefix"=prefixList, "pusedo_file"=res_files)
  return(res_df)
}

get_valid_path<-function(oldpath){
  result<-gsub('[ /:_()]+', '_', oldpath)
  return(result)
}

output_rawdata<-function(rawobj, outFile, Mtpattern, rRNApattern, hemoglobinPattern) {
  writeLines(rownames(rawobj), paste0(outFile, ".genes.txt"))

  saveRDS(rawobj, paste0(outFile, ".rawobj.rds"))

  png(paste0(outFile, ".top20.png"), width=3000, height=2000, res=300)
  par(mar = c(4, 8, 2, 1))
  C <- rawobj@assays$RNA@counts
  C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  mc<-rowMedians(C)
  most_expressed <- order(mc, decreasing = T)[20:1]
  tm<-as.matrix(Matrix::t(C[most_expressed,]))
  boxplot(tm, cex = 0.1, las = 1, xlab = "% total count per cell",
          col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
  dev.off()

  has_project = ifelse("project" %in% colnames(rawobj@meta.data), any(rawobj$orig.ident != rawobj$project),FALSE)
  has_sample = ifelse("sample" %in% colnames(rawobj@meta.data), any(rawobj$orig.ident != rawobj$sample), FALSE)

  draw_feature_qc(outFile, rawobj, "orig.ident")

  if(has_sample){
    draw_feature_qc(paste0(outFile, ".sample"), rawobj, "sample")
  }

  if(has_project){
    draw_feature_qc(paste0(outFile, ".project"), rawobj, "project")
  }

  rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
  rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]

  rawobj<-PercentageFeatureSet(object=rawobj, pattern=Mtpattern, col.name="percent.mt")
  rawobj<-PercentageFeatureSet(object=rawobj, pattern=rRNApattern, col.name = "percent.ribo")
  rawobj<-PercentageFeatureSet(object=rawobj, pattern=hemoglobinPattern, col.name="percent.hb")    

  draw_feature_qc(paste0(outFile, ".no_ribo"), rawobj, "orig.ident")

  if(has_sample){
    draw_feature_qc(paste0(outFile, ".no_ribo.sample"), rawobj, "sample")
  }

  if(has_project){
    draw_feature_qc(paste0(outFile, ".no_ribo.project"), rawobj, "project")
  }
}

XAxisRotation<-function(){
  return(theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}

save_session_info<-function(filename="sessionInfo.txt") {
  writeLines(capture.output(sessionInfo()), filename)
}

sub_cluster<-function(subobj, 
                     assay, 
                     by_sctransform, 
                     by_harmony, 
                     redo_harmony, 
                     curreduction, 
                     k_n_neighbors,
                     u_n_neighbors,
                     random.seed,
                     resolutions,
                     cur_npcs, 
                     cur_pca_dims,
                     vars.to.regress, 
                     essential_genes, 
                     key = "",
                     do_umap = TRUE,
                     reduction.name = "umap"
){
  n_half_cell=round(ncol(subobj) / 2)
  if(cur_npcs >= n_half_cell){
    cur_npcs = n_half_cell
    cur_pca_dims=1:cur_npcs
  }
  if(by_harmony){
    if(redo_harmony){
      if (length(unique(subobj$batch)) == 1){
        cat(key, "use old harmony result\n")
      }else{
        cat(key, "redo harmony\n")
        cat("RunPCA ... \n")
        subobj <- RunPCA(object = subobj, npcs=cur_npcs, assay=assay, verbose=FALSE)
        cat("RunHarmony ... \n")
        subobj <- RunHarmony(object = subobj,
                          assay.use = assay,
                          reduction = "pca",
                          dims.use = cur_pca_dims,
                          group.by.vars = "batch",
                          do_pca=FALSE)    
      }
    }else{
      #due to very limited cell numbers in small cluster, it may cause problem to redo harmony, 
      cat(key, "use old harmony result\n")
    }
    if(!("harmony" %in% names(subobj@reductions))){
      curreduction = "pca";
      cat(key, "no harmony, use PCA.\n")
    }
  }else{
    #https://github.com/satijalab/seurat/issues/5244
    if (by_sctransform) {
      cat(key, "use old sctransform result\n")
      #due to very limited cell numbers in small cluster, it may cause problem to redo sctransform at individual sample level, 
      #so we will keep the old data structure
      #subobj<-do_sctransform(subobj, vars.to.regress=vars.to.regress)
    }else{
      cat(key, "redo normalization\n")
      subobj<-do_normalization(subobj, selection.method="vst", nfeatures=3000, vars.to.regress=vars.to.regress, scale.all=FALSE, essential_genes=essential_genes)
    }
    cat(key, "RunPCA\n")
    subobj<-RunPCA(subobj, npcs=cur_npcs)
  }

  cat(key, "FindNeighbors\n")
  subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)

  cat(key, "FindClusters\n")
  subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolutions, verbose=FALSE)
  
  if(do_umap){
    cat(key, "RunUMAP\n")
    cur_min_dist = 0.3
    subobj<-RunUMAP(object = subobj, min.dist = cur_min_dist, reduction=curreduction, n.neighbors=u_n_neighbors, dims=cur_pca_dims, verbose = FALSE, reduction.name=reduction.name)
  }

  return(subobj)    
}

get_dot_height_vec<-function(vec){
  ngroups = length(unique(vec))
  result = max(1500, ngroups * 90 + 200)
  return(result)
}

get_dot_height<-function(obj, group.by){
  return(get_dot_height_vec(unlist(obj[[group.by]])))
}

output_barplot<-function(obj, sample_key, cell_key, filename){
  cts<-unlist(obj[[cell_key]])
  ct<-table(cts)
  ct<-ct[order(ct, decreasing = T)]
  ct_levels=names(ct)
  
  samples<-unlist(obj[[sample_key]])
  
  tbl<-data.frame(table(samples, cts))
  colnames(tbl)<-c("Sample", "Cell_type", "Cell_count")
  tbl$Cell_type<-factor(tbl$Cell_type, levels=ct_levels)
  
  g1<-ggplot(tbl, aes(Cell_type, Cell_count)) + geom_bar(aes(fill=Sample), stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("")
  g2<-ggplot(tbl, aes(Cell_type, Cell_count)) + geom_bar(aes(fill=Sample), position="fill", stat="identity") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  g<-g1+g2+plot_layout(ncol=1)
  png(filename, width=2000, height=3000, res=300)
  print(g)
  dev.off()
}

output_highlight_cell_by_cell_identity<-function(obj,
                                                 cell_identity,
                                                 prefix,
                                                 group.by="orig.ident",
                                                 name=group.by,
                                                 umap_width=2600,
                                                 cell_identity_umap="umap"){
  
  meta = obj@meta.data
  pcts = unlist(unique(obj[[cell_identity]]))
  for(pct in pcts){
    pct_str = celltype_to_filename(pct)
    cells<-rownames(meta)[meta[,cell_identity] == pct]
    subobj<-subset(obj, cells=cells)
    curprefix<-paste0(prefix, ".", cell_identity, ".", name, ".", pct_str)
    save_highlight_cell_plot(paste0(curprefix, ".cell.png"), subobj, group.by = group.by, reduction = cell_identity_umap, reorder=FALSE, title=pct)
  }
}

output_celltype_figures<-function(obj, 
                                  cell_identity, 
                                  prefix, 
                                  bubblemap_file, 
                                  cell_activity_database, 
                                  combined_ct_source, 
                                  group.by="orig.ident", 
                                  name=group.by, 
                                  umap_width=2600, 
                                  dot_width=4000, 
                                  cell_identity_order=NULL,
                                  all_umap="umap",
                                  cell_identity_umap="umap"){

  nct<-length(unique(unlist(obj[[cell_identity]])))

  if(is.null(cell_identity_order)){
    g<-get_dim_plot_labelby(obj = obj, 
                            label.by = cell_identity, 
                            title = "",
                            reduction = all_umap, 
                            legend.title = "")
  }else{
    g<-get_dim_plot(obj = obj,
                    group.by = cell_identity_order,
                    label.by = cell_identity, 
                    title = "",
                    reduction = all_umap, 
                    legend.title = "")
  }
  png(paste0(prefix, ".", cell_identity, ".umap.png"), width=umap_width, height=2000, res=300)
  print(g)
  dev.off()


  if(nct > 1){
    if(group.by != "batch"){
      if(is.null(cell_identity_order)){
        g<-get_bubble_plot( obj = obj, 
                            cur_res = NULL, 
                            cur_celltype = cell_identity, 
                            bubblemap_file = bubblemap_file, 
                            assay = "RNA", 
                            orderby_cluster = T)
      }else{
        g<-get_bubble_plot( obj = obj, 
                            cur_res = NULL, 
                            cur_celltype = cell_identity, 
                            bubblemap_file = bubblemap_file, 
                            assay = "RNA", 
                            orderby_cluster = FALSE,
                            group.by = cell_identity)
      }

      png(paste0(prefix, ".", cell_identity, ".dot.png"), width=dot_width, height=get_dot_height(obj, cell_identity), res=300)
      print(g)
      dev.off()
    }

    if(!all(is.na(cell_activity_database))){
      g<-get_celltype_marker_bubble_plot( obj = obj, 
                                          group.by = cell_identity, 
                                          cellType = cell_activity_database$cellType,
                                          weight = cell_activity_database$weight,
                                          n_markers = 5, 
                                          combined_ct_source=combined_ct_source)

      png(paste0(prefix, ".", cell_identity, ".ct_markers.bubbleplot.png"), width=dot_width, height=get_dot_height(obj, cell_identity), res=300)
      print(g)
      dev.off()
    }
  }

  output_barplot(obj, sample_key = group.by, cell_key = cell_identity, paste0(prefix, ".", cell_identity, ".", group.by, ".png"))

  save_highlight_cell_plot(filename = paste0(prefix, ".", cell_identity, ".cell.png"), 
                           obj = obj, 
                           group.by = cell_identity, 
                           reduction = all_umap,
                           reorder = is.null(cell_identity_order))

  output_highlight_cell_by_cell_identity(obj = obj,
                                         cell_identity = cell_identity,
                                         prefix = prefix,
                                         group.by = group.by,
                                         name = name,
                                         umap_width = umap_width,
                                         cell_identity_umap = cell_identity_umap)
  
  meta = obj@meta.data
  tb<-table(meta[,group.by], meta[,cell_identity])
  write.csv(tb, file=paste0(prefix, ".", cell_identity, ".", name,  "_celltype.csv"))

  mtb<-reshape2::melt(tb)
  colnames(mtb)<-c("Sample", "Celltype", "Cell")
  g1<-ggplot(mtb, aes(fill=Celltype, y=Cell, x=Sample)) + geom_bar(position="stack", stat="identity") + theme_bw3() + XAxisRotation() + xlab("")
  g2<-ggplot(mtb, aes(fill=Celltype, y=Cell, x=Sample)) + geom_bar(position="fill", stat="identity") + theme_bw3() + XAxisRotation() + xlab("")
  g<-g1+g2+plot_layout(ncol=1)
  png(paste0(prefix, ".", cell_identity, ".", name,  "_celltype.png"), width=3300, height=4000, res=300)
  print(g)
  dev.off()

  g1<-ggplot(mtb, aes(fill=Sample, y=Cell, x=Celltype)) + geom_bar(position="stack", stat="identity") + theme_bw3() + XAxisRotation() + xlab("")
  g2<-ggplot(mtb, aes(fill=Sample, y=Cell, x=Celltype)) + geom_bar(position="fill", stat="identity") + theme_bw3() + XAxisRotation() + xlab("")
  g<-g1+g2+plot_layout(ncol=1)
  png(paste0(prefix, ".", cell_identity, ".celltype_", name, ".png"), width=3300, height=4000, res=300)
  print(g)
  dev.off()
}

save_umap<-function(file_prefix, obj, umap_names=c("UMAP_1", "UMAP_2") ){
  umap<-FetchData(obj, umap_names)
  saveRDS(umap, paste0(file_prefix, ".rds"))
  write.csv(umap, paste0(file_prefix, ".csv"))
}

get_sig_gene_figure<-function(cell_obj, sigout, design_data, sig_gene, DE_by_cell=TRUE, is_between_cluster=FALSE, log_cpm=NULL){
  group_levels<-unique(design_data$Group)
  display_group_levels<-unique(design_data$DisplayGroup)

  groupColors<-c("blue", "red")
  names(groupColors)<-display_group_levels

  gmap<-unlist(split(design_data$Group, design_data$Sample))
  gdismap<-unlist(split(design_data$DisplayGroup, design_data$Sample))

  cell_obj$Group=factor(gmap[cell_obj$orig.ident], levels=group_levels)
  cell_obj$DisplayGroup=factor(gdismap[cell_obj$orig.ident], levels=display_group_levels)

  logFC<-sigout[sig_gene, "logFC"]
  FDR<-sigout[sig_gene,"FDR"]

  geneexp=FetchData(cell_obj,vars=c(sig_gene))
  colnames(geneexp)<-"Gene"
  colorRange<-c(min(geneexp), max(geneexp))
  fix.sc <- scale_color_gradientn(colors=c("lightgrey", "blue"), limits = colorRange)
  
  geneexp$Group<-cell_obj$DisplayGroup
  geneexp$Sample<-cell_obj$orig.ident
  
  title<-paste0(sig_gene, ' : logFC = ', round(logFC, 2), ", FDR = ", formatC(FDR, format = "e", digits = 2))
  
  if(is_between_cluster){
    p0<-ggplot(geneexp, aes(x=Group, y=Gene, col=Group)) + geom_violin() + geom_jitter(width = 0.2) + 
      facet_grid(~Sample) + theme_bw3() + 
      scale_color_manual(values = groupColors) +
      NoLegend() + xlab("") + ylab("Gene Expression")
    
    p1<-DimPlot(cell_obj, reduction = "umap", label=T, group.by="DisplayGroup") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5)) + xlim(xlim) + ylim(ylim)
    
    p2<-MyFeaturePlot(object = cell_obj, features=as.character(sig_gene), order=T)
    p<-p0+p1+p2+plot_layout(design="AA
BC")
    
  }else{
    p0<-ggplot(geneexp, aes(x=Sample, y=Gene, color=Group)) + 
      geom_violin() + geom_jitter(width = 0.2) + 
      theme_bw3() + theme_rotate_x_axis_label() +
      scale_color_manual(values = groupColors) +
      xlab("") + ylab("Gene Expression")
    
    if("subumap" %in% names(cell_obj@reductions)){
      p1<-MyFeaturePlot(object = cell_obj, features=sig_gene, order=T, reduction="subumap") + xlab("UMAP_1") + ylab("UMAP_2")
    }else{
      p1<-MyFeaturePlot(object = cell_obj, features=sig_gene, order=T, reduction="umap")
    }
    p1$data$Group=cell_obj@meta.data[rownames(p1$data), "DisplayGroup"]
    p1<-p1+facet_grid(~Group) + theme_bw3() + ggtitle("")
    
    if(!DE_by_cell){
      melt_cpm = as.data.frame(t(log_cpm[sig_gene,,drop=F]))
      colnames(melt_cpm)<-"CPM"
      melt_cpm$Group = factor(unlist(gdismap[rownames(melt_cpm)]), levels=display_group_levels)
      g0<-ggplot(melt_cpm, aes(x=Group, y=CPM, color=Group)) + geom_violin() + geom_boxplot(width=0.2) + geom_jitter(width = 0.1) +  
        theme_bw3() +
        scale_color_manual(values = groupColors) +
        xlab("") + ylab("CPM") + NoLegend()
      p<-p0+g0+p1+plot_layout(design="AAAAA
BCCCC")
    }else{
      p<-p0+p1+plot_layout(design="A
B")
    }
  }
  p<-p+ plot_annotation(title=title)
  return(p)
}

read_scrna_data<-function(fileName){
  if(dir.exists(fileName)){
    feature.names <- read.delim(paste0(fileName, "/features.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
    gene.column=ifelse(ncol(feature.names) > 1, 2, 1)
    counts = Read10X(fileName, gene.column=gene.column)
    if(is.list(counts)){
      if("protein_coding" %in% names(counts)){
        counts<-do.call(rbind, counts)
      }
    }
  } else if (grepl('.h5$', fileName)) {
    counts = Read10X_h5(fileName)
  } else if (grepl('.gz$', fileName)) {
    counts = data.frame(read_gzip_count_file(fileName, fileTitle, species))
  } else if (grepl('.rds$', fileName)) {
    counts = readRDS(fileName)
    if("Seurat" %in% class(counts)){
      counts = GetAssayData(counts, slot = "counts")
    }
  } else {
    stop(paste0("I don't know format of ", fileName))
  }
  
  adt.counts<-NULL
  if (is.list(counts) & ("Gene Expression" %in% names(counts))){
    adt.counts<-counts$`Antibody Capture`
    counts<-counts$`Gene Expression` 
  }
  return(list(counts=counts, adt.counts=adt.counts))
}

Plot_predictcelltype<-function(predict_celltype, topn=3, filename=NA, width=2000, height=2000) {
  library(heatmap3)

  cta_index<-apply(predict_celltype$cta,2,function(x){return(order(x,decreasing=T)[1:topn])})
  cta_index<-unique(sort(cta_index))
  cta_mat<- predict_celltype$cta[cta_index,]
  colnames(cta_mat)=paste0(colnames(cta_mat), ":", names(predict_celltype$max_cta))
  
  if(!is.na(filename)){
    png(filename, width=width, height=height, res=300)
  }
  heatmap3(cta_mat, scale="none", margin=c(10, 10), cexRow=1, cexCol=1, col=colorRampPalette(c("blue", "white", "red"))(n=1024))
  if(!is.na(filename)){
    dev.off()
  }
}

Plot_predictcelltype_ggplot2<-function(predict_celltype, topn=3, filename=NA, width=NA, height=NA, is_validation=FALSE) {
  library(ggplot2)
  library(scales)
  
  cta_index<-apply(predict_celltype$cta,2,function(x){return(order(x,decreasing=T)[1:topn])})
  cta_index<-unique(sort(cta_index))
  cta_mat<- predict_celltype$cta[cta_index,,drop=F]
  if(!is_validation){
    colnames(cta_mat)=paste0(colnames(cta_mat), ":", names(predict_celltype$max_cta))
  }

  n_ct=nrow(cta_mat)
  n_cluster=ncol(cta_mat)
  
  if(is.na(width)){
    width = max(1000, n_cluster * 50) + 500 
  }
  
  if(is.na(height)){
    height = max(1000, n_ct * 50) + 250 
  }
  
  ord <- hclust( dist(cta_mat, method = "euclidean"), method = "ward.D" )$order
  sorted_ct=rownames(cta_mat)[ord]

  ord <- hclust( dist(t(cta_mat), method = "euclidean"), method = "ward.D" )$order
  sorted_cluster=colnames(cta_mat)[ord]

  md = reshape2::melt(cta_mat)
  colnames(md)<-c("db_cell_type", "cluster", "cta_score")
  md$db_cell_type<-factor(md$db_cell_type, levels=sorted_ct)
  md$cluster<-factor(md$cluster, levels=sorted_cluster)

  colors = colorRampPalette(c("blue", "white", "red"))(n=1024)

  if(!is.na(filename)){
    png(filename, width=width, height=height, res=300)
  }
  g<-ggplot(md, aes(cluster, db_cell_type) ) +
    geom_tile(aes(fill = cta_score)) +
    scale_fill_gradientn(colours=colors) + theme_bw3(axis.x.rotate=TRUE) + xlab("") + ylab("")
  print(g)
  if(!is.na(filename)){
    dev.off()
  }
}
