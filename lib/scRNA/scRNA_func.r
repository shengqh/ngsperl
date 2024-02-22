library(harmony)
library(cowplot)
library(Seurat)
library(tools)
library(scales)
library(ggplot2)
library(patchwork)
library(Matrix.utils)
library(parallel)
library(data.table)
library(dplyr)
library(rlang)

#https://github.com/r-lib/lobstr/blob/main/R/mem.R
lobstr_node_size <- function() {
  bit <- 8L * .Machine$sizeof.pointer
  if (!(bit == 32L || bit == 64L)) {
    stop("Unknown architecture", call. = FALSE)
  }

  if (bit == 32L) 28L else 56L
}

lobstr_mem_used <- function() {
  rlang:::new_bytes(sum(gc()[, 1] * c(lobstr_node_size(), 8)))
}

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

get_value<-function(value, defaultValue){
  if(is.null(value)){
    return(defaultValue)
  }
  if(is.na(value)){
    return(defaultValue)
  }
  if(value == ""){
    return(defaultValue)
  }
  return(value)
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
    stop(paste0("file not exists: ", filepath))
  }
  
  if(file.info(filepath)$size == 0){
    return(TRUE)
  }

  return(FALSE)
}

is_seurat_object<-function(obj){
  return("Seurat" %in% class(obj))
}

read_file_map<-function(file_list_path, sep="\t", header=F, do_unlist=TRUE){
  if(grepl('.csv$', file_list_path)){
    tbl<-read.csv(file_list_path, header=header)
  }else{
    tbl<-read.table(file_list_path, sep=sep, header=header)
  }
  result<-split(tbl$V1, tbl$V2)
  if(do_unlist){
    result<-unlist(result)
  }
  return(result)
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

theme_rotate_x_axis_label <- function(angle=90, vjust=0.5, hjust=1) {
  theme(axis.text.x = element_text(angle=angle, vjust=vjust, hjust=hjust))
}

theme_bw3 <- function (axis.x.rotate=F, angle=90, vjust=0.5, hjust=1) { 
  is_ggplot2_newver = packageVersion("ggplot2") >= "3.4.0"

  if(is_ggplot2_newver){
    eline = element_line(colour = "black", linewidth = 0.5)
  }else{
    eline = element_line(colour = "black", size = 0.5)
  }

	result = theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, colour = 'black'),
      panel.border = element_rect(fill = NA, color = "black"),			
      axis.line = eline
    )
  if (axis.x.rotate){
    result = result + theme_rotate_x_axis_label(angle=angle,vjust=vjust,hjust=hjust)
  }

  return(result)
}

get_heatmap_height<-function(ngenes){
  result<-max(3000, min(20000, ngenes * 60 + 1000))
  return(result)
}

get_heatmap_width<-function(nclusters){
  result<-max(3000, min(10000, nclusters * 300 + 1000))
  return(result)
}

MyDimPlot<-function(...){
  #set raster=FALSE 
  #When number of cells larger than 100000, it might throw error: 
  #Problem while converting geom to grob, Error occurred in the 1st layer.
  g<-DimPlot(raster=FALSE, repel=TRUE, ...) + theme(aspect.ratio=1)
  return(g)
}

#https://github.com/satijalab/seurat/issues/1836
#For visualization, using sctransform data is also fine.
MyDoHeatMap<-function(obj, max_cell=5000, ...){
  if(ncol(obj) > max_cell){
    subsampled <- subset(obj, cells = sample(colnames(obj), size=max_cell, replace=F))
    g<-DoHeatmap(subsampled, ...)
  }else{
    g<-DoHeatmap(obj, ...)
  }
  return(g)
}

MyFeaturePlot<-function(object, assay="RNA", ...){
  old_assay=DefaultAssay(object)
  DefaultAssay(object)=assay
  g=FeaturePlot(object=object, ...) + theme(aspect.ratio=1)
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

do_normalization<-function(obj, selection.method="vst", nfeatures=2000, vars.to.regress=NULL, scale.all=FALSE, essential_genes=NULL) {
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

do_sctransform<-function(rawobj, vars.to.regress=NULL, return.only.var.genes=FALSE, mc.cores=1, use_sctransform_v2=FALSE) {
  vst.flavor = ifelse(use_sctransform_v2, "v2", "v1")

  print(paste0("performing SCTransform by ", vst.flavor, " ..."))
  nsamples=length(unique(rawobj$orig.ident))
  if(nsamples > 1){
    mc.cores = check_mc_cores(mc.cores)

    print("  split objects ...")
    objs<-SplitObject(object = rawobj, split.by = "orig.ident")
    rm(rawobj)

    print("  perform sctransform ...")
    if(mc.cores > 1){
      objs<-mclapply(objs, function(x){
        print(paste0("    sctransform ", unique(x$orig.ident), " ..."))
        x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE, vst.flavor=vst.flavor)
        return(x)
      }, mc.cores=mc.cores)  
    }else{
      objs<-lapply(objs, function(x){
        print(paste0("    sctransform ", unique(x$orig.ident), " ..."))
        x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE, vst.flavor=vst.flavor)
        return(x)
      })  
    }
    print("  sctransform done")

    print("  merge samples ...")
    obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
    #https://github.com/satijalab/seurat/issues/2814
    VariableFeatures(obj[["SCT"]]) <- rownames(obj[["SCT"]]@scale.data)
    rm(objs)
    return(obj)
  }else{
    print("  perform sctransform ...")
    rawobj<-SCTransform(rawobj, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE, vst.flavor=vst.flavor)
    print("  sctransform done")
    return(rawobj)
  }
}

do_harmony<-function(obj, by_sctransform, vars.to.regress, has_batch_file, batch_file, pca_dims, essential_genes=NULL, mc.cores=1,use_sctransform_v2=TRUE){
  if(by_sctransform){
    #now perform sctranform
    obj<-do_sctransform(obj, vars.to.regress=vars.to.regress,mc.cores=mc.cores,use_sctransform_v2=use_sctransform_v2)
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

init_cutoffs<-function(all_samples, myoptions, filter_config_file=""){
  ##cutoff dataframe for each sample to filter empty droplets
  if(filter_config_file != ""){
    Cutoffs<-fread(filter_config_file, data.table=F, header=TRUE)
    missed_samples = all_samples[!(all_samples %in% Cutoffs$sample)]
    if(length(missed_samples) > 0){
      missed_Cutoffs<-data.frame( sample=missed_samples,
                                  nFeature_cutoff_min=myoptions$nFeature_cutoff_min ,
                                  nFeature_cutoff_max=myoptions$nFeature_cutoff_max,
                                  nCount_cutoff=myoptions$nCount_cutoff, 
                                  mt_cutoff=myoptions$mt_cutoff, 
                                  cluster_remove=c(""),
                                  stringsAsFactors = F)
      Cutoffs<-rbind(Cutoffs, missed_Cutoffs)
    }
    Cutoffs = Cutoffs[Cutoffs$sample %in% all_samples,,drop=FALSE]
  }else{
    Cutoffs<-data.frame(
      sample=all_samples,
      nFeature_cutoff_min=myoptions$nFeature_cutoff_min ,
      nFeature_cutoff_max=myoptions$nFeature_cutoff_max,
      nCount_cutoff=myoptions$nCount_cutoff, 
      mt_cutoff=myoptions$mt_cutoff, 
      cluster_remove=c(""),
      stringsAsFactors = F)
  }
  rownames(Cutoffs)<-Cutoffs$sample
  Cutoffs$nFeature_cutoff_min<-as.numeric(Cutoffs$nFeature_cutoff_min)
  Cutoffs$nFeature_cutoff_max<-as.numeric(Cutoffs$nFeature_cutoff_max)
  Cutoffs$nCount_cutoff<-as.numeric(Cutoffs$nCount_cutoff)
  Cutoffs$mt_cutoff<-as.numeric(Cutoffs$mt_cutoff)
  return(Cutoffs)
}

preprocessing_rawobj<-function(rawobj, myoptions, prefix, filter_config_file=""){
  Mtpattern= myoptions$Mtpattern
  rRNApattern=myoptions$rRNApattern
  
  Cutoffs=init_cutoffs(
    all_samples=unique(rawobj$orig.ident), 
    myoptions=myoptions, 
    filter_config_file=filter_config_file)

  rawCells<-data.frame(table(rawobj$orig.ident))
  
  plot1 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
    geom_hline(data=Cutoffs, aes(yintercept=mt_cutoff, color=sample))  + 
    geom_vline(data=Cutoffs, aes(xintercept=nCount_cutoff, color=sample)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10))

  plot2 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_hline(data=Cutoffs, aes(yintercept=nFeature_cutoff_min, color=sample)) + 
    geom_hline(data=Cutoffs, aes(yintercept=nFeature_cutoff_max, color=sample)) + 
    geom_vline(data=Cutoffs, aes(xintercept=nCount_cutoff, color=sample)) 

  p<-plot1+plot2
  ggsave(paste0(prefix, ".qc.1.png"), p, width=11, height=5, dpi=300, units="in", bg="white")

  mt<-data.frame(mt=rawobj$percent.mt, Sample=rawobj$orig.ident, nFeature=log10(rawobj$nFeature_RNA), nCount=log10(rawobj$nCount_RNA))

  nsample<-length(unique(mt$Sample))
  ncol=ceiling(sqrt(nsample))
  nrow=ceiling(nsample/ncol)
  width=min(10000, ncol * 1200)
  height=min(10000, nrow * 1000)

  vcutoffs = Cutoffs
  vcutoffs$log10_nCount_cutoff = log10(vcutoffs$nCount_cutoff)
  vcutoffs$log10_nFeature_cutoff_min = log10(vcutoffs$nFeature_cutoff_min)
  vcutoffs$log10_nFeature_cutoff_max = log10(vcutoffs$nFeature_cutoff_max)
  vcutoffs$Sample = vcutoffs$sample

  g1<-ggplot(mt, aes(y=mt,x=nCount) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_hline(data=vcutoffs, aes(yintercept=mt_cutoff, color=Sample)) + 
    geom_vline(data=vcutoffs, aes(xintercept=log10_nCount_cutoff, color=Sample)) +
    ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    facet_wrap(~Sample, ncol=ncol, nrow=nrow) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  
  ggsave(paste0(prefix, ".qc.2.png"), g1, width=width, height=height, dpi=300, units="px", bg="white")
  
  g1<-ggplot(mt, aes(y=mt,x=nFeature) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    geom_hline(data=vcutoffs, aes(yintercept=mt_cutoff, color=Sample)) + 
    geom_vline(data=vcutoffs, aes(xintercept=log10_nFeature_cutoff_min, color=Sample)) +
    geom_vline(data=vcutoffs, aes(xintercept=log10_nFeature_cutoff_max, color=Sample)) +
    ylab("Percentage of mitochondrial") + xlab("log10(number of feature)") +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    facet_wrap(~Sample, ncol=ncol, nrow=nrow) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))

  ggsave(paste0(prefix, ".qc.3.png"), g1, width=width, height=height, dpi=300, units="px", bg="white")

  finalList<-list()

  #filter cells
  meta=rawobj@meta.data
  cells = c()
  for(idx in c(1:nrow(Cutoffs))){
    sub_meta=meta[meta$orig.ident==Cutoffs$sample[idx],]
    filteredsub_meta<-subset(sub_meta, nFeature_RNA >= Cutoffs$nFeature_cutoff_min[idx] & 
                                  nFeature_RNA <= Cutoffs$nFeature_cutoff_max[idx] & 
                                  nCount_RNA >= Cutoffs$nCount_cutoff[idx] & 
                                  percent.mt <= Cutoffs$mt_cutoff[idx])
    cells<-c(cells, rownames(filteredsub_meta))
  }
  rawobj<-subset(rawobj, cells=cells)

  filteredCells<-data.frame(table(rawobj$orig.ident))
  qcsummary<-merge(rawCells, filteredCells, by = "Var1")
  colnames(qcsummary)<-c("Sample", "RawCell", "ValidCell")
  qcsummary$DiscardCell<-qcsummary$RawCell-qcsummary$ValidCell
  qcsummary$DiscardRate<-qcsummary$DiscardCell / qcsummary$RawCell
  qcsummary=merge(Cutoffs, qcsummary, by.x="sample", by.y="Sample")
  write.csv(qcsummary, file=paste0(prefix, ".filtered.cell.csv"), row.names=F)
  
  g<-VlnPlot(object = rawobj, features = c("percent.mt", "nFeature_RNA", "nCount_RNA"), group.by="orig.ident")
  ggsave(paste0(prefix, ".qc.4.png"), g, width=12, height=5, dpi=300, units="in", bg="white")
  
  finalList$filter<-qcsummary
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

  if("ADT" %in% names(obj)){
    adt_names=rownames(obj$ADT@counts)
    writeLines(adt_names, paste0(outFile, ".ADT.txt"))
    for(adt in adt_names){
      g<-FeaturePlot(obj, features=paste0("adt_", adt), cols = c("lightgrey", "red"), min.cutoff=0.01, max.cutoff=0.99, order=TRUE) + ggtitle(paste0(adt, " protein")) + theme(aspect.ratio=1)
      ggsave(paste0(outFile, ".", adt, ".png"), g, width=5, height=5, dpi=300, units="in", bg="white")
    }

    common_genes=intersect(adt_names, rownames(obj@assays$RNA@counts))
    for(adt in common_genes){
      g1<-FeaturePlot(obj, features=paste0("adt_", adt), cols = c("lightgrey", "red"), min.cutoff=0.01, max.cutoff=0.99, order=TRUE) + ggtitle(paste0(adt, " protein")) + theme(aspect.ratio=1)
      g2<-FeaturePlot(obj, features=paste0("rna_", adt), cols = c("lightgrey", "red"), min.cutoff=0.01, max.cutoff=0.99, order=TRUE) + ggtitle(paste0(adt, " RNA")) + theme(aspect.ratio=1)
      g<-g1+g2+plot_layout(ncol=2)
      ggsave(paste0(outFile, ".", adt, ".common.png"), g, width=10, height=5, dpi=300, units="in", bg="white")
    }
  }
}

do_read_bubble_genes<-function(bubblemap_file, allgenes=c(), species="Hs"){
  library("readxl")
  library("tidyr")

  if(grepl(".txt$", bubblemap_file)){
    genes <- read.table(bubblemap_file, header=F, sep="\t", stringsAsFactors = F)
    if(any(grepl(",", genes$V2))){
      genes <- genes[,c(2,1)]
    }
  }else{
    genes <- data.frame(read_xlsx(bubblemap_file, sheet = 1))
    if(colnames(genes)[1] %in% c('Count', 'Index')){
      genes <- genes[,c(2:ncol(genes))]
    }

    while(all(is.na(genes[,1]))){
      genes <- genes[,c(2:ncol(genes))]
    }

    if(colnames(genes)[2] == 'Marker.Gene'){
      genes <- genes[,c(2:ncol(genes))]
    }
  }

  colnames(genes)[1:2] = c("gene", "cell_type")  

  #don't use the marker genes other than thhe first column, it might just annotation
  #if you want to use those genes, put into the first column
  # if(ncol(genes) > 2){
  #   #for some excel files, the third column is addtional genes
  #   if(any(grepl(',', genes[,1])) & any(grepl(',', genes[,3]))){
  #     genes[,3][is.na(genes[,3])]<-""
  #     genes[,1] = paste0(genes[,1], ",", genes[,3])
  #   }

  #   genes<-genes[,c(1,2)]
  # }

  for(idx in c(2:nrow(genes))){
    if(is.na(genes[idx,"cell_type"])){
      genes[idx,"cell_type"]=genes[idx-1,"cell_type"]
    }
  }

  #for some excel file, multiple genes are in same row but seperated by ',', 
  #use sepearte_rows to put them in different rows.
  genes<-data.frame(separate_rows(genes, gene, sep="[, ]+"))
  genes<-genes[genes$gene != "",]

  if(!is.null(species)){
    if(tolower(species) == "mm" | tolower(species) == "mouse"){
      genes$gene = toMouseGeneSymbol(genes$gene)
    }
  }

  gene_names=genes$gene
  gene_names[gene_names=="PECAM"] = "PECAM1"
  gene_names[gene_names=="HGD1B"] = "HGD"
  gene_names[gene_names=="EpCAM"] = "EPCAM"
  gene_names[gene_names=="CD25"] = "IL2RA"
  gene_names[gene_names=="ACTAA2"] = "ACTA2"
  gene_names[gene_names=="MTND6"] = "MT-ND6"
  gene_names[gene_names=="FOXJ!"] = "FOXJ1"
  
  genes$gene<-gene_names

  genes<-genes[!duplicated(genes),]

  if(length(allgenes) > 0){
    miss_genes=setdiff(genes$gene, allgenes)
    writeLines(miss_genes, con="miss_gene.csv")

    genes<-genes[genes$gene %in% allgenes,]
  }else{
  }
  genes$cell_type=factor(genes$cell_type, levels=unique(genes$cell_type))
  
  return(genes)
}

read_bubble_genes<-function(bubblemap_files, allgenes=c(), species="Hs"){
  result = NULL
  bubblemap_file=bubblemap_files[1]
  for(bubblemap_file in bubblemap_files){
    genes = do_read_bubble_genes(bubblemap_file, allgenes, species)
    if(is.null(result)){
      result = genes
    }else{
      result = rbind(result, genes)
    }
  }
  result$cell_type=factor(result$cell_type, levels=unique(result$cell_type))
  
  return(result)
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

get_dot_plot<-function(obj, group.by, gene_groups, assay="RNA", rotate.title=TRUE, use_blue_yellow_red=TRUE ){
  genes=unique(unlist(gene_groups))
  assaydata=GetAssayData(obj, assay=assay, slot="data")
  if(!all(genes %in% rownames(assaydata))){
    missed_genes = genes[!(genes %in% rownames(assaydata))]
    missed_genes=missed_genes[c(1:min(5, length(missed_genes)))]
    stop(paste0("some genes are not in ", assay, " assay, here is the first few:", paste0(missed_genes, collapse = ",")))
  }

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
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = "Identity") +
    theme_cowplot() + 
    facet_grid(facets = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") + 
    theme(panel.spacing = unit(x = 1,units = "lines"), strip.background = element_blank())
    
  if(use_blue_yellow_red){
    plot <- plot + scale_colour_gradient2(low="blue", mid="yellow", high="red", midpoint=0 )
  }else{
    plot <- plot + scale_color_gradient(low="lightgray", high="blue")
  }
  
  g=plot + xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                             axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                                             strip.background = element_blank())
  
  if(rotate.title){
    g=g+theme(strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
  }
  return(g)
}

get_dot_width<-function(g, min_width=4400){
  if(!all(c("features.plot","feature.groups") %in% colnames(g$data))){
    stop(paste0("features.plot or feature.groups is not in ", paste0(colnames(g$data), collapse = ",")))
  }
  ngenes = nrow(g$data[!duplicated(g$data[,c("features.plot","feature.groups")]),])
  ngroups = length(unique(g$data$feature.groups))
  width=ngenes * 50 + ngroups * 40 + 400
  return(max(width, min_width))
}

get_dot_height_num<-function(ngroups, min_height=2000, height_per_entry=80, height_additional_space=1000){
  result = max(min_height, ngroups * height_per_entry + height_additional_space)
  return(result)
}

get_dot_height_vec<-function(vec){
  ngroups = length(unique(vec))
  return(get_dot_height_num(ngroups))
}

get_dot_height<-function(obj, group.by){
  return(get_dot_height_vec(unlist(obj[[group.by]])))
}

get_bubble_plot<-function(obj, 
                          cur_res, 
                          cur_celltype, 
                          bubblemap_file, 
                          assay="RNA", 
                          orderby_cluster=FALSE, 
                          split.by=NULL, 
                          rotate.title=TRUE, 
                          group.by=NULL, 
                          use_blue_yellow_red=TRUE, 
                          species="Hs"){
                            
  old_assay = DefaultAssay(obj)
  DefaultAssay(obj) = assay
  allgenes=rownames(obj)
  DefaultAssay(obj) = old_assay
  
  if(assay=="ADT"){
    genes_df <- read_bubble_genes(bubblemap_file, allgenes, species=NULL)
    if(nrow(genes_df) == 0){
      stop(paste0("no proteins in ", bubblemap_file," found in data, please double check your protein names. They should match with following protein names: ", paste0(allgenes, collapse = ", ")))
    }
  }else{
    genes_df <- read_bubble_genes(bubblemap_file, allgenes, species=species)
    if(nrow(genes_df) == 0){
      stop(paste0("no genes in ", bubblemap_file," found in data, is it posible wrong species used in function? species=", species))
    }
  }

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

  g<-get_dot_plot(obj, group.by, gene_groups, assay, rotate.title=rotate.title, use_blue_yellow_red=use_blue_yellow_red)
  
  return(g)
}

get_sub_bubble_plot<-function(obj, obj_res, subobj, subobj_res, bubblemap_file, add_num_cell=FALSE, species=NULL, assay="RNA"){
  old_meta<-obj@meta.data
  
  obj$fake_layer=paste0("fake_", unlist(obj@meta.data[,obj_res]))

  sr = as.character(subobj@meta.data[,subobj_res])

  if(add_num_cell){
    num_tbl = table(sr)
    num_tbl = num_tbl[order(num_tbl, decreasing = T)]
    sub_levels = paste0(names(num_tbl), " (", num_tbl, ")")
    obj@meta.data[colnames(subobj), "fake_layer"] = paste0(sr, " (", num_tbl[sr], ")")
  }else{
    obj@meta.data[colnames(subobj), "fake_layer"] = sr
    sub_levels = levels(subobj@meta.data[, subobj_res])
  }

  g<-get_bubble_plot( obj, 
                      assay=assay,
                      group.by="fake_layer", 
                      bubblemap_file = bubblemap_file, 
                      species=species)

  obj@meta.data=old_meta
  
  g$data<-g$data[!grepl("^fake_", g$data$id),]
  if(all(!is.null(sub_levels))){
    g$data$id<-factor(g$data$id, levels=sub_levels)
  }

  return(g)
}

draw_bubble_plot<-function(obj, cur_res, cur_celltype, bubble_map_file, prefix, width=5500, height=3000, rotate.title=TRUE, species="Hs"){
  g<-get_bubble_plot( obj, 
                      cur_res, 
                      cur_celltype, 
                      bubble_map_file, 
                      rotate.title=rotate.title, 
                      species=species)
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
  
read_object<-function(obj_file, meta_rds=NULL, columns=NULL, sample_name=NULL){
  obj=readRDS(obj_file)
  if(is.list(obj)){
    if(!is.null(sample_name)){
      if(sample_name %in% names(obj)){
        obj<-obj[[sample_name]]
      }
    }
  }

  if(is.list(obj)){
    if("obj" %in% names(obj)){
      obj=obj$obj
    }else{
      stop(paste0("No obj found in list with names: ", paste0(names(obj), collapse=", ")))
    }
  }
  
  if(!is_seurat_object(obj)){
    stop(paste0("It is not Seurat object: ", obj_file))
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
        if(column %in% colnames(meta.data)){
          obj=AddMetaData(obj, meta.data[,column], col.name=column)
        }
      }
    }
  }
  return(obj)
}

read_object_from_file_list<-function(file_list_path, meta_rds=NULL, columns=NULL){
  df=fread(file_list_path, header=F)
  sample_name=df$V2[1]
  obj_file=df$V1[1]
  obj=read_object(obj_file, meta_rds=meta_rds, columns=columns, sample_name=sample_name)
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

  if(!("project" %in% colnames(rawobj@meta.data))){
    rawobj$project = rawobj$orig.ident
  }

  if(!("sample" %in% colnames(rawobj@meta.data))){
    rawobj$sample = rawobj$orig.ident
  }

  if(ident_name == "orig.ident"){
    if(any(rawobj$orig.ident != rawobj$sample)){
      meta<-rawobj@meta.data
      meta$os = paste0(meta$orig.ident, ":", meta$sample)
      meta<-meta[!duplicated(meta$os),,drop=F]
      smap=split(as.character(meta$sample), as.character(meta$orig.ident))
      smap2=lapply(smap,function(x){
        paste0(x, collapse=",")
      })
      ct$Source=unlist(smap2[ct$Sample])
    }

    if(any(rawobj$project != rawobj$sample)){
      meta<-rawobj@meta.data
      meta$os = paste0(meta$orig.ident, ":", meta$project)
      meta<-meta[!duplicated(meta$os),,drop=F]
      smap=split(as.character(meta$project), as.character(meta$orig.ident))
      smap2=lapply(smap,function(x){
        paste0(x, collapse=",")
      })
      ct$Project=unlist(smap2[ct$Sample])
    }
  }
  write.table(ct, paste0(prefix, ".cell.txt"), sep="\t", row.names=F)
  
  if("Project" %in% colnames(ct)){
    g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Project))
  }else{
    g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Sample)) + NoLegend()
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
  ldata<-as.character(unlist(obj[[layer]]))
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

  g<-MyDimPlot(obj, group.by=group.by, label=label, reduction=reduction, split.by=split.by, ...)+ 
    scale_color_manual(legend.title, values=scolors, labels = cts, guide = guide_legend(ncol=ncol)) + 
    ggtitle(title)
  return(g)
}

build_dummy_cluster<-function(obj, label.by, new_cluster_name, new_cluster_name_label=paste0(new_cluster_name, "_label")){
  groups<-as.character(obj@meta.data[,label.by])
  gt<-table(groups)
  gt<-gt[order(gt, decreasing=T)]
  dummy_cluster<-c(0:(length(gt)-1))
  names(dummy_cluster)<-names(gt)
  dc<-factor(dummy_cluster[groups], levels=dummy_cluster)
  obj@meta.data[,new_cluster_name]<-dc
  obj@meta.data[,new_cluster_name_label]<-paste0(obj@meta.data[,new_cluster_name], ": ", groups)
  return(obj)
}

get_dim_plot_labelby<-function(obj, label.by, title=label.by, label=T, legend.title=label.by, reduction="umap", split.by=NULL, ncol=1, label_has_cluster=FALSE, ...){
  group.by="dummy_cluster"
  group.label="dummy_label"

  obj<-build_dummy_cluster(obj, label.by, group.by, group.label)

  group.label=ifelse(label_has_cluster, label.by, group.label)

  g<-get_dim_plot(obj, group.by=group.by, label.by=group.label, label=label, title=title, legend.title=legend.title, reduction=reduction, split.by=split.by, ncol=ncol, ...)
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
    g0<-MyDimPlot(obj, label=F, cells.highlight =cells, reduction = reduction) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", ct)) + NoLegend()
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
  if("seurat_clusters" %in% colnames(clusterDf)){
    cts = as.character(unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = T), cluster_name]))
  }else{
    tbl = table(clusterDf[,cluster_name])
    tbl = tbl[order(tbl, decreasing = T)]
    cts = names(tbl)
  }
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
  if("ADT" %in% names(rawobj)){
    writeLines(rownames(rawobj$ADT@counts), paste0(outFile, ".antibodies.txt"))
  }

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

  # rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
  # rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]

  # rawobj<-PercentageFeatureSet(object=rawobj, pattern=Mtpattern, col.name="percent.mt", assay="RNA")
  # rawobj<-PercentageFeatureSet(object=rawobj, pattern=rRNApattern, col.name = "percent.ribo", assay="RNA")
  # rawobj<-PercentageFeatureSet(object=rawobj, pattern=hemoglobinPattern, col.name="percent.hb", assay="RNA")    

  # draw_feature_qc(paste0(outFile, ".no_ribo"), rawobj, "orig.ident")

  # if(has_sample){
  #   draw_feature_qc(paste0(outFile, ".no_ribo.sample"), rawobj, "sample")
  # }

  # if(has_project){
  #   draw_feature_qc(paste0(outFile, ".no_ribo.project"), rawobj, "project")
  # }
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

  if(by_harmony & !("harmony" %in% names(subobj@reductions))){
    redo_harmony=TRUE
  }

  if(redo_harmony){
    curreduction = "harmony";
    if(!("batch" %in% colnames(subobj))){
      subobj$batch = subobj$orig.ident
    }

    if (length(unique(subobj$batch)) == 1){
      cat(key, "use old harmony result\n")
    }else{
      cat(key, "redo harmony\n")
      subobj <- RunHarmony(object = subobj,
                        assay.use = assay,
                        dims.use = cur_pca_dims,
                        group.by.vars = "batch",
                        do_pca=TRUE)    
    }
    if(!("harmony" %in% names(subobj@reductions))){
      curreduction = "pca";
      cat(key, "no harmony, use PCA.\n")
    }
  }else if (by_harmony) {
    cat(key, "use old harmony result\n")
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
  cat("curreduction =", curreduction, "\n")

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
                                  cell_identity_umap="umap",
                                  species="Hs"){

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
  ggsave(paste0(prefix, ".", cell_identity, ".umap.png"), g, width=umap_width, height=2000, dpi=300, units="px", bg="white")
  rm(g)

  if(nct > 1){
    if(group.by != "batch"){
      if(is.null(cell_identity_order)){
        g<-get_bubble_plot( obj = obj, 
                            cur_res = NULL, 
                            cur_celltype = cell_identity, 
                            bubblemap_file = bubblemap_file, 
                            assay = "RNA", 
                            orderby_cluster = T,
                            species=species)
      }else{
        g<-get_bubble_plot( obj = obj, 
                            cur_res = NULL, 
                            cur_celltype = cell_identity, 
                            bubblemap_file = bubblemap_file, 
                            assay = "RNA", 
                            orderby_cluster = FALSE,
                            group.by = cell_identity,
                            species=species)
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

  if(!is_between_cluster){
    ddata<-design_data[!duplicated(design_data$Sample),]

    gmap<-unlist(split(ddata$Group, ddata$Sample))
    gdismap<-unlist(split(ddata$DisplayGroup, ddata$Sample))

    cell_obj@meta.data$Group=factor(gmap[cell_obj$orig.ident], levels=group_levels)
    cell_obj@meta.data$DisplayGroup=factor(gdismap[cell_obj$orig.ident], levels=display_group_levels)
  }

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
    p0<-ggplot(geneexp, aes(x=Group, y=Gene, col=Group)) + geom_violin() + geom_jitter(width = 0.2)

    if(length(unique(geneexp$Sample)) > 1){
      p0 = p0 + facet_grid(~Sample)
    }
    p0 = p0 + theme_bw3() + 
      scale_color_manual(values = groupColors) +
      NoLegend() + xlab("") + ylab("Gene Expression")
    
    p1<-MyDimPlot(cell_obj, reduction = "umap", label=T, group.by="DisplayGroup") + NoLegend() + ggtitle("Cluster") + theme(plot.title = element_text(hjust=0.5)) + xlim(xlim) + ylim(ylim)
    
    p2<-MyFeaturePlot(object = cell_obj, features=as.character(sig_gene), order=T, raster=FALSE)
    p<-p0+p1+p2+plot_layout(design="AA
BC")
    
  }else{
    p0<-ggplot(geneexp, aes(x=Sample, y=Gene, color=Group)) + 
      geom_violin() + geom_jitter(width = 0.2) + 
      theme_bw3() + theme_rotate_x_axis_label() +
      scale_color_manual(values = groupColors) +
      xlab("") + ylab("Gene Expression")
    
    if("subumap" %in% names(cell_obj@reductions)){
      p1<-MyFeaturePlot(object = cell_obj, features=sig_gene, order=T, reduction="subumap", raster=FALSE) + xlab("UMAP_1") + ylab("UMAP_2")
    }else{
      p1<-MyFeaturePlot(object = cell_obj, features=sig_gene, order=T, reduction="umap", raster=FALSE)
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

read_scrna_data<-function(fileName, keep_seurat=FALSE){
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
    counts = tryCatch ({
        Read10X_h5(fileName)
      },
      error = function(e) {
        message(paste("Read file as 10X_h5 failed, try as CellBender format:", fileName))
        library(scCustomize)
        Read_CellBender_h5_Mat(file_name=fileName)
      }
    )
  } else if (grepl('.gz$', fileName)) {
    counts = data.frame(read_gzip_count_file(fileName, fileTitle, species))
  } else if (grepl('.rds$', fileName)) {
    counts = readRDS(fileName)
    if(is_seurat_object(counts)){
      if(!keep_seurat){
        counts = GetAssayData(counts, slot = "counts")
      }
    }
  } else {
    stop(paste0("I don't know format of ", fileName))
  }
  
  adt.counts<-NULL
  if (is.list(counts) & ("Gene Expression" %in% names(counts))){
    adt.counts<-counts$`Antibody Capture`
    counts<-counts$`Gene Expression` 
  }
  return(list(counts=ceiling(counts), adt.counts=adt.counts))
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
    width = max(1500, n_cluster * 60) + 1000 
  }
  
  if(is.na(height)){
    height = max(1000, n_ct * 60) + 400 
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

match_cell_names<-function(sample_name, source_meta, target_meta){
  cur_cells = intersect(rownames(source_meta), rownames(target_meta))

  if(length(cur_cells) > 0){
    return(source_meta)
  }

  cur_names<-paste0(sample_name, "_", rownames(source_meta))
  cur_cells = intersect(cur_names, rownames(target_meta))

  if(length(cur_cells) > 0){
    rownames(source_meta)<-cur_names
    return(source_meta)
  }

  cur_names<-paste0(sample_name, "_", rownames(target_meta))
  cur_cells = intersect(rownames(source_meta), cur_names)
  if(length(cur_cells) > 0){
    rownames(source_meta)<-gsub(paste0(sample_name, "_"), "", rownames(source_meta))
    return(source_meta)
  }

  if("project" %in% colnames(target_meta)){
    obj_meta = target_meta[target_meta$project == sample_name,]
  }else{
    obj_meta = target_meta[target_meta$sample == sample_name,]
  }

  if(!"orig.cell" %in% colnames(target_meta)){
    obj_meta$orig.cell<-gsub(".+_", "", rownames(obj_meta))
  }

  cur_cells = intersect(rownames(source_meta), obj_meta$orig.cell)
  if(length(cur_cells) == 0){
    stop(paste0("I don't know how to map source meta cell ", rownames(source_meta)[1], " with target meta ", rownames(target_meta)[1]))
  }
  
  matched = rownames(source_meta) %in% obj_meta$orig.cell
  old_names = rownames(source_meta)[matched]
  cell_map = unlist(split(rownames(obj_meta), obj_meta$orig.cell))
  rownames(source_meta)[matched] = unlist(cell_map[old_names])

  return(source_meta)
}

fill_meta_info<-function(sample_name, source_meta, target_meta, source_columns, target_column, is_character=TRUE){
  source_meta<-match_cell_names(sample_name, source_meta, target_meta)
  cur_cells = intersect(rownames(source_meta), rownames(target_meta))
  for(source_column in source_columns){
    if(source_column %in% colnames(source_meta)){
      if(is_character){
        target_meta[cur_cells, target_column] = as.character(source_meta[cur_cells, source_column])
      }else{
        target_meta[cur_cells, target_column] = source_meta[cur_cells, source_column]
      }
      break
    }
  }
  return(target_meta)
}

fill_meta_info_list<-function(source_meta_file_list, target_meta, source_columns, target_column, is_character=TRUE){
  source_map<-read_file_map(source_meta_file_list)
  cur_name=names(source_map)[1]
  for(cur_name in names(source_map)){
    source_meta_file=source_map[cur_name]
    source_meta=readRDS(source_meta_file)
    target_meta=fill_meta_info(cur_name, source_meta, target_meta, source_columns, target_column, is_character)
  }
  return(target_meta)
}

# The default FeaturePlot function in Seurat doesn't handle the order correctly. We need to fix it.
my_feature_plot<-function(obj, gene, high_color="red", umap1 = "UMAP_1", umap2 = "UMAP_2", split.by=NULL, point.size=0.2, order=TRUE){
  if(!is.null(split.by)){
    gdata<-FetchData(obj, c(umap1, umap2, gene, split.by))
  }else{
    gdata<-FetchData(obj, c(umap1, umap2, gene))
  }

  if(order){
    gdata<-gdata[order(gdata[,3]),]
  }

  if(max(gdata[,3]) == 0){
    high_color = "lightgray"
  }

  g<-ggplot(gdata, aes(!!sym(umap1), !!sym(umap2), color=!!sym(gene))) + 
    geom_point(size=point.size) + 
    scale_color_gradient(low="lightgray", high=high_color) + 
    theme_bw3() +
    theme(aspect.ratio=1)

  if(!is.null(split.by)){
    g<-g + facet_grid(as.formula(paste0("~ ", split.by))) +
      theme(strip.background = element_blank())
  }

  g
}

get_barplot<-function(
  ct_meta, 
  bar_file=NULL, 
  cluster_name="display_layer", 
  validation_columns=c("orig.ident","SignacX","SingleR"), 
  calc_height_per_cluster=200, 
  calc_width_per_cell=50){

  valid_columns = intersect(validation_columns, colnames(ct_meta))
  if(length(valid_columns) == 0){
    stop(paste0("No column found in meta: ", paste0(validation_columns, collapse=", ")))
  }
  
  alltbl=NULL
  col_name="SignacX"
  for(col_name in valid_columns){
    tbl = data.frame(table(ct_meta[,cluster_name], ct_meta[,col_name]))
    v1 = as.numeric(as.character(tbl$Var1))
    if(all(is.na(v1))){
      v1 = as.character(tbl$Var1)
    }
    tbl$Var1 = v1
    tbl$Category=col_name

    alltbl<-rbind(alltbl, tbl)
  }

  g<-ggplot(alltbl, aes(Var2, Freq, fill=Var2)) + 
    geom_bar(width=0.5, stat = "identity") + 
    facet_grid(Var1~Category, scales = "free", space='free_x') + 
    theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend() +
    theme(strip.text.y = element_text(angle = 0))

  if(!is.null(bar_file)){
    height = max(1000, length(unique(alltbl$Var1)) * calc_height_per_cluster + 500)
    width = max(1000, length(unique(alltbl$Var2)) * calc_width_per_cell) + 400

    ggsave(bar_file, g, width=width, height=height, dpi=300, units="px", bg="white")
  }
  
  return(g)
}

h5ad_to_h5seurat <- function(h5ad_file){
  library(Seurat)
  library(SeuratData)
  library(SeuratDisk)
  library(hdf5r)

  Convert(h5ad_file, dest = "h5seurat")
  h5seurat_file <- gsub(".h5ad$", ".h5seurat", h5ad_file)

  #https://github.com/mojaveazure/seurat-disk/issues/109
  f <- H5File$new(h5seurat_file, "r+")
  groups <- f$ls(recursive = TRUE)

  for (name in groups$name[grepl("categories", groups$name)]) {
    names <- strsplit(name, "/")[[1]]
    names <- c(names[1:length(names) - 1], "levels")
    new_name <- paste(names, collapse = "/")
    f[[new_name]] <- f[[name]]
  }

  for (name in groups$name[grepl("codes", groups$name)]) {
    names <- strsplit(name, "/")[[1]]
    names <- c(names[1:length(names) - 1], "values")
    new_name <- paste(names, collapse = "/")
    f[[new_name]] <- f[[name]]
    grp <- f[[new_name]]
    grp$write(args = list(1:grp$dims), value = grp$read() + 1)
  }

  f$close_all()
}

read_object_from_rawfile<-function(sample_name, file_path, species, ensembl_map=NULL){
  cat("reading", sample_name, ":", file_path, "\n")
  lst = read_scrna_data(file_path, keep_seurat=TRUE)

  counts = lst$counts  
  adt.counts = lst$adt.counts

  if(is_seurat_object(counts)){
    sobj<-counts
    rm(counts)

    sobj$orig.ident = sample_name
  }else{
    rs<-rowSums(counts)
    counts<-counts[rs>0,]

    if(!is.null(ensembl_map)){
      gtf_counts<-counts[!(rownames(counts) %in% names(ensembl_map)),]
      ensembl_counts<-counts[(rownames(counts) %in% names(ensembl_map)),]
      gene_names<-unlist(ensembl_map[rownames(ensembl_counts)])
      gene_counts<-DelayedArray::rowsum(ensembl_counts, gene_names)
      counts<-rbind(gtf_counts, gene_counts)
    }

    if (species=="Mm") {
      rownames(counts)<-toMouseGeneSymbol(rownames(counts))
    }
    if (species=="Hs") {
      rownames(counts)<-toupper(rownames(counts))
    }
    sobj = CreateSeuratObject(counts = counts, project = sample_name)
    sobj$orig.ident <- sample_name
    rm(counts)
  }
  return(sobj)
}


iterate_celltype<-function(obj, 
                           previous_celltypes, 
                           previous_layer, 
                           iter_name, 
                           cur_layermap, 
                           npcs, 
                           resolution, 
                           random.seed, 
                           by_sctransform, 
                           by_harmony, 
                           curprefix, 
                           iter, 
                           vars.to.regress,
                           bubblemap_file, 
                           essential_genes,
                           species="Hs"){
  meta = obj@meta.data
  
  assay=ifelse(by_sctransform, "SCT", "RNA")

  files = get_empty_files()
  
  all_cur_cts<-NULL
  pct<-previous_celltypes[length(previous_celltypes)]
  
  #previous_celltypes<-c("Platelets")
  for(pct in previous_celltypes){
    pct_str = celltype_to_filename(pct)

    key = paste0("iter", iter, ": ", pct, ":")
    cells<-rownames(meta)[meta[,previous_layer] == pct]
    if(length(cells) == 0){#no cell left for this cell type
      next
    }
    
    subobj<-subset(obj, cells=cells)

    subobj[["oumap"]] = subobj[["umap"]]
    
    stopifnot(all(subobj[[previous_layer]] == pct))
    
    pca_npcs<-min(round(length(cells)/2), 50)
    cur_npcs=min(pca_npcs, npcs)
    cur_pca_dims=1:cur_npcs

    k_n_neighbors<-min(cur_npcs, 20)
    u_n_neighbors<-min(cur_npcs, 30)

    DefaultAssay(subobj)<-assay

    curreduction=ifelse(by_harmony, "harmony", "pca")

    if(pct != "Unassigned") {
      subobj = sub_cluster(subobj = subobj, 
                            assay =  assay, 
                            by_sctransform = by_sctransform, 
                            by_harmony = by_harmony, 
                            redo_harmony = redo_harmony,
                            curreduction = curreduction, 
                            k_n_neighbors = k_n_neighbors,
                            u_n_neighbors = u_n_neighbors,
                            random.seed = random.seed,
                            resolutions = resolution,
                            cur_npcs = cur_npcs, 
                            cur_pca_dims = cur_pca_dims,
                            vars.to.regress = vars.to.regress,
                            essential_genes = essential_genes
                            )
    }else{
      cat(key, "FindNeighbors\n")
      subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)

      cat(key, "FindClusters\n")
      subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolution, verbose=FALSE)
    }
    
    cat(key, "Cell type annotation\n")
    cur_cts<-subobj[[previous_layer]]

    cluster<-"seurat_clusters"
    data_norm=get_seurat_average_expression(subobj, cluster)
    
    predict_celltype<-ORA_celltype(data_norm,cell_activity_database$cellType,cell_activity_database$weight)
    
    cta_rds_file=paste0(curprefix, ".", pct_str, ".cta.rds")
    saveRDS(predict_celltype, cta_rds_file)
    files<-rbind(files, c(previous_layer, iter_name, pct, "cta_rds", cta_rds_file))

    if(length(predict_celltype$max_cta) > 1){
      cta_png_file=paste0(curprefix, ".", pct_str, ".cta.png")
      Plot_predictcelltype_ggplot2( predict_celltype, 
                            filename=cta_png_file)
      files<-rbind(files, c(previous_layer, iter_name, pct, "cta_png", cta_png_file))
    }

    new_cluster_ids<-names(predict_celltype$max_cta)
    
    cur_cts<-cbind_celltype(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts)
    stopifnot(all(colnames(subobj) == rownames(cur_cts)))
    subobj = AddMetaData(subobj, cur_cts$seurat_clusters, "seurat_clusters")
    subobj = AddMetaData(subobj, cur_cts$cell_type, "cell_type")
    subobj = AddMetaData(subobj, cur_cts$seurat_cell_type, "seurat_cell_type")
    subobj = AddMetaData(subobj, cur_cts$raw_cell_type, "raw_cell_type")
    subobj = AddMetaData(subobj, cur_cts$raw_seurat_cell_type, "raw_seurat_cell_type")
    
    #using RNA assay for visualization
    DefaultAssay(subobj)<-assay

    if(pct == "Unassigned"){
      g0<-MyDimPlot(obj, label=F, group.by=previous_layer) + ggtitle(pct)
    }else{
      g0<-MyDimPlot(obj, label=F, cells.highlight =cells) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))
    }
    g1<-get_dim_plot_labelby(subobj, reduction="oumap", label.by="cell_type") + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("New cell type in old UMAP")
    g2<-get_dim_plot(subobj, reduction="oumap", group.by="seurat_clusters", label.by="raw_seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat raw cell type in old UMAP")
    g3<-get_dim_plot(subobj, reduction="oumap", group.by="seurat_clusters", label.by="seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat cell type in old UMAP")
    g<-g0+g1+g2+g3+plot_layout(nrow=2)
    umap_celltype_file = paste0(curprefix, ".", pct_str, ".old_umap.png")
    png(umap_celltype_file, width=4600, height=4000, res=300)
    print(g)
    dev.off()

    files<-rbind(files, c(previous_layer, iter_name, pct, "old_umap", umap_celltype_file))

    if(pct != "Unassigned"){
      g1<-get_dim_plot(subobj, group.by="seurat_clusters", label.by="raw_seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Raw cell type in new UMAP")
      g2<-get_dim_plot(subobj, group.by="seurat_clusters", label.by="seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat cell type in new UMAP")

      g<-g1+g2+plot_layout(nrow=1)
      umap_cluster_file = paste0(curprefix, ".", pct_str, ".new_umap.png")
      png(umap_cluster_file, width=4600, height=2000, res=300)
      print(g)
      dev.off()
      files<-rbind(files, c(previous_layer, iter_name, pct, "new_umap", umap_cluster_file))
    }

    if(pct == "Unassigned"){
      g<-get_bubble_plot( obj = subobj, 
                          bubblemap_file = bubblemap_file, 
                          group.by = "raw_seurat_cell_type",
                          species=species)
    }else{
      g<-get_sub_bubble_plot( obj, 
                              previous_layer, 
                              subobj, 
                              "raw_seurat_cell_type", 
                              bubblemap_file,
                              species=species)
    }
    dot_file = paste0(curprefix, ".", pct_str, ".dot.png")
    png(dot_file, width=get_dot_width(g), height=get_dot_height(subobj, cluster), res=300)
    print(g)
    dev.off()
    files<-rbind(files, c(previous_layer, iter_name, pct, "dot", dot_file))

    all_cur_cts<-rbind(all_cur_cts, cur_cts)
    
    # if(previous_layer == "layer0"){
    #   obj[['umap']] = subobj[['umap']]
    # }
    rm(subobj)
  }
  colnames(files)<-colnames(get_empty_files())
  return(list("all_cur_cts"=all_cur_cts, "files"=files))
}

layer_cluster_celltype<-function(obj, 
                                 previous_layer, 
                                 cur_layer, 
                                 cur_layermap, 
                                 npcs, 
                                 resolution, 
                                 random.seed, 
                                 by_sctransform, 
                                 by_harmony, 
                                 prefix, 
                                 vars.to.regress,
                                 bubblemap_file,
                                 essential_genes,
                                 species=species){
  meta<-obj@meta.data
  
  previous_celltypes<-unique(meta[[previous_layer]])

  files = get_empty_files()

  if(length(previous_celltypes) > 0){
    cat("cluster and annotate cell types:", paste0(previous_celltypes, ", "))
    iter = 1
    while(TRUE){
      cat("Iteration ", iter, "\n")
      
      iter_name=paste0("iter", iter)

      previous_celltypes<-previous_celltypes[order(previous_celltypes)]
      
      curprefix = paste0(prefix, ".iter", iter)

      iter_meta_file = paste0(curprefix, ".csv")
      iter_meta_rds = paste0(curprefix, ".rds")

      cat("  Call iterate_celltype ...\n")
      lst<-iterate_celltype(obj, 
                            previous_celltypes, 
                            previous_layer, 
                            iter_name, 
                            cur_layermap, 
                            npcs, 
                            resolution, 
                            random.seed, 
                            by_sctransform, 
                            by_harmony, 
                            curprefix, 
                            iter, 
                            vars.to.regress,
                            bubblemap_file, 
                            essential_genes,
                            species=species)

      all_cur_cts<-lst$all_cur_cts
      cur_files<-lst$files

      files<-rbind(files, cur_files)
      
      stopifnot(all(rownames(all_cur_cts) %in% colnames(obj)))
      
      iter_clusters = paste0(iter_name, "_clusters")
      iter_raw = paste0(iter_name, "_raw")

      obj[[iter_name]] = obj[[previous_layer]]
      obj[[iter_clusters]] = obj[[paste0(previous_layer, "_clusters")]]
      obj[[iter_raw]] = obj[[paste0(previous_layer, "_raw")]]

      obj@meta.data[rownames(all_cur_cts), iter_name]=unlist(all_cur_cts[, "cell_type"])
      obj@meta.data[rownames(all_cur_cts), iter_clusters]=unlist(as.character(all_cur_cts[, "seurat_clusters"]))
      obj@meta.data[rownames(all_cur_cts), iter_raw]=unlist(all_cur_cts[, "raw_cell_type"])
      
      pre_disagree<-all_cur_cts[all_cur_cts[, previous_layer] != all_cur_cts[,'cell_type'],,drop=F]
      if(nrow(pre_disagree) > 0){
        cat("Found unmatched cell type\n")
        print(table(pre_disagree[,previous_layer], pre_disagree$cell_type))
        write.csv(obj@meta.data, iter_meta_file)
        saveRDS(obj@meta.data, iter_meta_rds)
        previous_celltypes = unique(c(unlist(pre_disagree[,previous_layer]), unlist(pre_disagree$cell_type)))
        previous_celltypes<-previous_celltypes[previous_celltypes %in% unlist(unique(obj[[iter_name]]))]
        previous_layer = iter_name
        iter = iter + 1
        next
      }else{
        obj[[cur_layer]] = obj[[iter_name]]
        obj[[paste0(cur_layer, "_clusters")]] = obj[[iter_clusters]]
        obj[[paste0(cur_layer, "_raw")]] = obj[[iter_raw]]
        write.csv(obj@meta.data, iter_meta_file)
        saveRDS(obj@meta.data, iter_meta_rds)
        break
      }
    }
  }
  
  #using RNA assay for visualization
  DefaultAssay(obj)<-"RNA"

  g<-get_dim_plot_labelby(obj, label.by = cur_layer, label=T)

  png(paste0(prefix, ".", cur_layer, ".umap.png"), width=3300, height=3000, res=300)
  print(g)
  dev.off()

  if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
    g2<-get_bubble_plot(obj, 
                        NULL, 
                        cur_layer, 
                        bubblemap_file, 
                        assay="RNA",
                        species=species)
    png(paste0(prefix, ".", cur_layer, ".dot.png"), width=get_dot_width(g2), height=get_dot_height(obj, cur_layer), res=300)
    print(g2)
    dev.off()
  }

  write.csv(obj@meta.data, file=paste0(prefix, ".", cur_layer, ".meta.csv"))
  saveRDS(obj@meta.data, file=paste0(prefix, ".", cur_layer, ".meta.rds"))

  write.csv(all_cur_cts, file=paste0(prefix, ".", cur_layer, ".details.csv"))
  
  return(list(obj=obj, files=files))
}


do_analysis<-function(tmp_folder,
                      cur_folder,
                      obj, 
                      layer2map, 
                      npcs, 
                      resolution, 
                      random.seed, 
                      by_sctransform, 
                      by_harmony, 
                      prefix, 
                      vars.to.regress, 
                      bubblemap_file, 
                      essential_genes,
                      by_individual_sample,
                      species="Hs" ) {
  setwd(tmp_folder)
  reslist1<-layer_cluster_celltype( obj = obj,
                                    previous_layer = "layer0", 
                                    cur_layer = "layer4", 
                                    cur_layermap = layer2map, 
                                    npcs = npcs, 
                                    resolution = resolution, 
                                    random.seed = random.seed, 
                                    by_sctransform = by_sctransform, 
                                    by_harmony = by_harmony, 
                                    prefix = prefix, 
                                    vars.to.regress = vars.to.regress,
                                    bubblemap_file = bubblemap_file,
                                    essential_genes = essential_genes,
                                    species=species)
  obj=reslist1$obj
  files=reslist1$files
  rm(reslist1)

  setwd(cur_folder)

  write.csv(files, paste0(prefix, ".iter_png.csv"))

  celltypes<-unique(obj$layer4)
  celltypes<-celltypes[order(celltypes)]
  ctdf<-data.frame("celltype"=celltypes, "resolution"=0.01)
  write.table(ctdf, paste0(prefix, ".scDynamic.celltype_res.txt"), row.names=F, sep="\t", quote=F)

  obj<-factorize_layer(obj, "layer4")
  Idents(obj)<-"layer4"

  saveRDS(obj@meta.data, paste0(prefix, ".scDynamic.meta.rds"))
  write.csv(obj@meta.data, paste0(prefix, ".scDynamic.meta.csv"))

  save_umap(paste0(prefix, ".scDynamic.umap"), obj)

  if(length(unique(obj$layer4)) > 1){
    #find markers for all cell types
    all_markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    all_top10<-get_top10_markers(all_markers)
    all_top10<-unique(all_top10$gene)

    obj<-myScaleData(obj, all_top10, "RNA")
    g<-MyDoHeatMap(obj, max_cell=5000, assay="RNA", features = all_top10, group.by = "layer4", angle = 90) + NoLegend()

    width<-get_heatmap_width(length(unique(Idents(obj))))
    height<-get_heatmap_height(length(all_top10))
    png(paste0(prefix, ".layer4.heatmap.png"), width=width, height=height, res=300)
    print(g)
    dev.off()
  }
  
  output_celltype_figures(obj, 
    "layer4", 
    prefix, 
    bubblemap_file, 
    cell_activity_database, 
    combined_ct_source, 
    group.by="orig.ident", 
    name="sample",
    species=species)

  if(!by_individual_sample){
    #output individual sample dot plot, with global scaled average gene expression.
    obj$sample_layer4<-paste0(obj$orig.ident, ":", obj$layer4)
    g<-get_bubble_plot( obj = obj, 
                        cur_res = NA, 
                        cur_celltype = "sample_layer4", 
                        bubblemap_file = bubblemap_file, 
                        assay = "RNA", 
                        orderby_cluster = F,
                        species=species)
    gdata<-g$data
    gdata$id<-gsub(".+: ","",gdata$id)
    gdata$sample<-gsub(":.+","",gdata$id)
    gdata$id<-gsub(".+:","",gdata$id)

    for(sample in unique(gdata$sample)){
      sdata<-gdata[gdata$sample == sample,]
      g$data=sdata
      dot_file = paste0(prefix, ".layer4.", sample, ".dot.png")
      png(dot_file, width=get_dot_width(g), height=get_dot_height(obj, "layer4"), res=300)
      print(g)
      dev.off()
    }

    has_batch<-FALSE
    if("batch" %in% colnames(obj@meta.data)){
      if("sample" %in% colnames(obj@meta.data)){
        has_batch=any(obj$batch != obj$sample)
      }else{
        has_batch=any(obj$batch != obj$orig.ident)
      }
    }
    if(has_batch){
      output_celltype_figures(obj, "layer4", prefix, bubblemap_file, cell_activity_database, combined_ct_source, group.by="batch", name="batch")
    }
  }

  obj$seurat_layer4=paste0(obj$layer4_clusters, ": ", obj$layer4_raw)

  cur_celltype="layer4"
  for(pct in unique(unlist(obj[[cur_celltype]]))){
    cells=colnames(obj)[obj[[cur_celltype]] == pct]
    subobj=subset(obj, cells=cells)
    subobj$seurat_layer4=paste0(subobj$layer4_clusters, ": ", subobj$layer4_raw)
    g<-get_dim_plot(subobj, group.by="layer4_clusters", label.by="seurat_layer4", reduction="umap", legend.title="")

    png(paste0(prefix, ".", cur_celltype, ".", celltype_to_filename(pct), ".umap.png"), width=2400, height=2000, res=300)
    print(g)
    dev.off()
  }

  if(by_individual_sample){
    tb=data.frame("cell"=colnames(obj), "cell_type"=obj$layer4, category=prefix)
    output_file=paste0(prefix,".dynamic.html")
    rmdfile = "seurat_scDynamic_one_layer_one_resolution.rmd"
    rmarkdown::render(rmdfile, output_file=output_file)
    return(list(html=output_file, ct_count=tb))
  }else{
    return(obj)
  }
}

factor_by_count<-function(vec){
  tbl=table(vec)
  tbl=tbl[tbl > 0]
  tbl=tbl[order(tbl, decreasing=T)]
  res=factor(vec, levels=names(tbl))
  return(res)
}
