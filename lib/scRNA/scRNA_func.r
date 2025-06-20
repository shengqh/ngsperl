
load_install<-function(library_name, library_sources=library_name){
  if(!require(library_name, character.only = T)){
    BiocManager::install(library_sources, ask=FALSE)
  }
  library(library_name, character.only = T)
}

load_install("harmony")
load_install("cowplot")
load_install("Seurat")
load_install("tools")
load_install("scales")
load_install("ggplot2")
load_install("patchwork")
load_install("Matrix.utils", "cvarrichio/Matrix.utils")
load_install("parallel")
load_install("data.table")
load_install("plyr")
load_install("dplyr")
load_install("rlang")
load_install("scCustomize")
load_install("SeuratData", "satijalab/seurat-data")
load_install("SeuratWrappers", "satijalab/seurat-wrappers")
load_install("BiocParallel")

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

ggvenn_2_set<-function(venn_data){
  library(ggvenn)
  g=ggvenn(venn_data, set_name_size=0) + 
    annotate("text", x = -0.8, y = 1.15, label = names(venn_data)[1], size = 5, hjust = 0.5) +
    annotate("text", x = 0.8, y = 1.15, label = names(venn_data)[2], size = 5, hjust = 0.5) 
  return(g)
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
  if(hasArg(raster)){
    g<-DimPlot(...) + theme(aspect.ratio=1)
  }else{
    g<-DimPlot(raster=FALSE, ...) + theme(aspect.ratio=1)
  }
  return(g)
}

#https://github.com/satijalab/seurat/issues/1836
#For visualization, using sctransform data is also fine.
MyDoHeatMap<-function(obj, max_cell=5000, features=NULL, assay=NULL, slot = "scale.data", ...){
  if(ncol(obj) > max_cell){
    cur_obj <- subset(obj, cells = sample(colnames(obj), size=max_cell, replace=F))
  }else{
    cur_obj <- obj
  }
  if(slot=="scale.data"){
    scale_data=MyGetAssayData(cur_obj, assay=assay, slot="scale.data")
    if(!all(features %in% rownames(scale_data))){
      cur_obj <- ScaleData(cur_obj, features=features)
    }
  }
  g<-DoHeatmap(cur_obj, features=features, assay=assay, slot=slot, ...)
  return(g)
}

MyFeaturePlot<-function(object, assay="RNA", use_scCustom=TRUE, ...){
  old_assay=DefaultAssay(object)
  DefaultAssay(object)=assay
  if(use_scCustom){
    g=FeaturePlot_scCustom(object, ...) + theme(aspect.ratio=1)
  }else{
    g=FeaturePlot(object, ...) + theme(aspect.ratio=1)
  }
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
  counts=MyGetAssayData(curobj,assay="RNA",slot="counts")
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

is_seurat_5_plus<-function(obj){
  return(Version(obj) > '5')
}

#sometimes, only the assay will be converted to Seurat 5, so we need to check the assay
is_assay_5_plus<-function(obj, assay){
  return(class(obj[[assay]]) == 'Assay5')
}

has_data<-function(obj, assay, slot){
  if(is_assay_5_plus(obj, assay)){
    return(slot %in% names(obj@assays[[assay]]@layers))
  }else{
    return(slot %in% names(obj@assays[[assay]]))
  }
}

MyGetAssayData<-function(obj, assay, slot){
  if(is_seurat_5_plus(obj)){
    return(GetAssayData(obj, assay=assay, layer=slot))
  # if(is_assay_5_plus(obj, assay)){
  #   cur_assay=obj[[assay]]
  #   return(LayerData(cur_assay, layer=slot))
  }else{
    return(GetAssayData(obj, assay=assay, slot=slot))
  }
}

do_normalization<-function( obj, 
                            selection.method="vst", 
                            nfeatures=2000, 
                            vars.to.regress=NULL, 
                            scale.all=FALSE, 
                            essential_genes=NULL,
                            ignore_variable_genes=NULL) {
  DefaultAssay(obj)<-"RNA"

  cat("NormalizeData ... \n")
  obj <- NormalizeData(obj, verbose = FALSE)
  cat("NormalizeData done ... \n")
  
  cat("FindVariableFeatures ... \n")
  obj <- FindVariableFeatures(obj, selection.method = selection.method, nfeatures = nfeatures, verbose = FALSE)
  cat("FindVariableFeatures done ... \n")

  if(length(ignore_variable_genes) > 0){
    cat("Remove ignore variable genes ... \n")
    vgenes=VariableFeatures(obj)
    VariableFeatures(obj) <- setdiff(vgenes, ignore_variable_genes)
    cat("Remove some variable genes done ... \n")
  }
  
  if(scale.all){
    features = rownames(obj)
  }else{
    features=VariableFeatures(obj)
    if (length(essential_genes) > 0){
      features = unique(c(features, essential_genes))
    }
  }
  cat("ScaleData: ", length(features), "genes ...\n")
  obj <- ScaleData(obj, vars.to.regress=vars.to.regress, features=features, verbose = FALSE)
  cat("ScaleData done ... \n")

  stopifnot(dim(MyGetAssayData(obj, assay="RNA", slot="scale.data"))[1] > 0)

  return(obj)
}

do_sctransform<-function( rawobj, 
                          vars.to.regress=NULL, 
                          return.only.var.genes=FALSE, 
                          mc.cores=1, 
                          use_sctransform_v2=TRUE,
                          ignore_variable_genes=NULL) {
  vst.flavor = ifelse(use_sctransform_v2, "v2", "v1")

  print(paste0("performing SCTransform by ", vst.flavor, " ..."))
  nsamples=length(unique(rawobj$orig.ident))
  if(nsamples > 1){
    mc.cores = check_mc_cores(mc.cores)

    print("  SplitObject ...")
    objs<-SplitObject(object = rawobj, split.by = "orig.ident")
    rm(rawobj)

    print("  perform sctransform ...")
    if(mc.cores > 1){
      objs<-mclapply(objs, function(x){
        print(paste0("    SCTransform ", unique(x$orig.ident), " ..."))
        x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE, vst.flavor=vst.flavor)
        return(x)
      }, mc.cores=mc.cores)  
    }else{
      objs<-lapply(objs, function(x){
        print(paste0("    SCTransform ", unique(x$orig.ident), " ..."))
        x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE, vst.flavor=vst.flavor)
        return(x)
      })  
    }
    print("  sctransform done")

    print("  merge samples ...")
    #collapse: If ‘TRUE’, merge layers of the same name together; if ‘FALSE’, appends ‘labels’ to the layer name
    #collapse layers is not supported yet
    obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")

    if(is_assay_5_plus(obj, "RNA")){
      print("  JoinLayers ...")
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    }

    #https://github.com/satijalab/seurat/issues/2814
    sct_var_genes=rownames(obj[["SCT"]]@scale.data)
    VariableFeatures(obj[["SCT"]]) <- setdiff(sct_var_genes, ignore_variable_genes)
    rm(objs)
    return(obj)
  }else{
    print("  perform sctransform ...")
    rawobj<-SCTransform(rawobj, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=return.only.var.genes, verbose = FALSE, vst.flavor=vst.flavor)
    print("  sctransform done")
    sct_var_genes=rownames(obj[["SCT"]]@scale.data)
    VariableFeatures(rawobj[["SCT"]]) <- setdiff(sct_var_genes, ignore_variable_genes)
    return(rawobj)
  }
}

do_harmony<-function(obj, by_sctransform, 
                    vars.to.regress, 
                    has_batch_file, 
                    batch_file, 
                    pca_dims, 
                    essential_genes=NULL, 
                    mc.cores=1,
                    use_sctransform_v2=TRUE,
                    ignore_variable_genes=NULL) {
  if(by_sctransform){
    #now perform sctranform
    obj<-do_sctransform(obj, 
                        vars.to.regress=vars.to.regress,
                        mc.cores=mc.cores,
                        use_sctransform_v2=use_sctransform_v2,
                        ignore_variable_genes=ignore_variable_genes)
    assay="SCT"
  }else{
    assay="RNA"
  }

  #no matter if we will use sctransform, we need normalized RNA assay for visualization and cell type annotation
  #data slot for featureplot, dotplot, cell type annotation and scale.data slot for heatmap
  #we need to do sctransform first, then do RNA assay normalization and scale, otherwise, after sctransform, the scale.data slot will be discarded.
  obj<-do_normalization(obj, 
                        selection.method="vst", 
                        nfeatures=2000, 
                        vars.to.regress=vars.to.regress, 
                        scale.all=FALSE, 
                        essential_genes=essential_genes,
                        ignore_variable_genes=ignore_variable_genes)

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

preprocessing_rawobj<-function(rawobj, myoptions, prefix, filter_config_file="", assay="RNA"){
  Mtpattern= myoptions$Mtpattern
  rRNApattern=myoptions$rRNApattern
  hemoglobinPattern=myoptions$hemoglobinPattern
  
  all_samples=unique(rawobj$orig.ident)
  Cutoffs=init_cutoffs(
    all_samples=all_samples, 
    myoptions=myoptions, 
    filter_config_file=filter_config_file)

  Remove_hemoglobin=is_one(myoptions$Remove_hemoglobin)
  Remove_rRNA=is_one(myoptions$Remove_rRNA)
  Remove_MtRNA=is_one(myoptions$Remove_MtRNA)

  remove_genes = c()
  if(Remove_MtRNA){
    mt_genes <- rownames(rawobj)[grepl(Mtpattern, rownames(rawobj))]
    cat("Remove mt genes: ", length(mt_genes), "\n")
    remove_genes = c(remove_genes, mt_genes)
  }
  if(Remove_rRNA){
    rRNA_genes <- rownames(rawobj)[grepl(rRNApattern, rownames(rawobj))]
    cat("Remove rRNA genes: ", length(rRNA_genes), "\n")
    remove_genes = c(remove_genes, rRNA_genes)
  }
  if(Remove_hemoglobin && !is.null(hemoglobinPattern) && hemoglobinPattern != ""){
    hemoglobin_genes <- rownames(rawobj)[grepl(hemoglobinPattern, rownames(rawobj))]
    cat("Remove hemoglobin genes: ", length(hemoglobin_genes), "\n")
    remove_genes = c(remove_genes, hemoglobin_genes)
  }
  if(length(remove_genes) > 0){
    cat("Remove total genes: ", length(remove_genes), "\n")
    new_genes = setdiff(rownames(rawobj), remove_genes)
    rawobj = subset(rawobj, features = new_genes)
  }
  rawobj<-PercentageFeatureSet(object=rawobj, pattern=Mtpattern, col.name="percent.mt", assay=assay)
  rawobj<-PercentageFeatureSet(object=rawobj, pattern=rRNApattern, col.name = "percent.ribo", assay=assay)
  if(!is.null(hemoglobinPattern) && hemoglobinPattern != ""){
    rawobj<-PercentageFeatureSet(object=rawobj, pattern=hemoglobinPattern, col.name="percent.hb", assay=assay)    
  }

  rawCells<-data.frame(table(rawobj$orig.ident))
  
  readCol = paste0("nCount_", assay)
  featureCol = paste0("nFeature_", assay)

  plot1 <- FeatureScatter(object = rawobj, feature1 = readCol, feature2 = "percent.mt") + 
    geom_hline(data=Cutoffs, aes(yintercept=mt_cutoff, color=sample))  + 
    geom_vline(data=Cutoffs, aes(xintercept=nCount_cutoff, color=sample)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10))

  plot2 <- FeatureScatter(object = rawobj, feature1 = readCol, feature2 = featureCol) +
    geom_hline(data=Cutoffs, aes(yintercept=nFeature_cutoff_min, color=sample)) + 
    geom_hline(data=Cutoffs, aes(yintercept=nFeature_cutoff_max, color=sample)) + 
    geom_vline(data=Cutoffs, aes(xintercept=nCount_cutoff, color=sample)) 

  p<-plot1+plot2

  if(length(all_samples) > 10){
    p<-p & NoLegend()
  }
  ggsave(paste0(prefix, ".qc.1.png"), p, width=11, height=5, dpi=300, units="in", bg="white")

  mt<-data.frame( mt=rawobj$percent.mt, 
                  Sample=rawobj$orig.ident, 
                  nFeature=log10(rawobj@meta.data[, featureCol]), 
                  nCount=log10(rawobj@meta.data[, readCol]))

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
    filteredsub_meta<-subset(sub_meta, sub_meta[,featureCol] >= Cutoffs$nFeature_cutoff_min[idx] & 
                                  sub_meta[,featureCol] <= Cutoffs$nFeature_cutoff_max[idx] & 
                                  sub_meta[,readCol] >= Cutoffs$nCount_cutoff[idx] & 
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
  
  g<-VlnPlot(object = rawobj, features = c("percent.mt", featureCol, readCol), group.by="orig.ident")
  ggsave(paste0(prefix, ".qc.4.png"), g, width=12, height=5, dpi=300, units="in", bg="white")
  
  finalList$filter<-qcsummary
  finalList$rawobj<-rawobj
  return(finalList)
}

output_integration_dimplot<-function(obj, outFile, has_batch_file, qc_genes=NULL){
  g<-FeaturePlot_scCustom(obj, features="percent.mt") + ggtitle("Percentage of mitochondrial genes")
  width=1700
  ncol=1
  
  if("percent.hb" %in% colnames(obj@meta.data)){
    g2<-FeaturePlot_scCustom(obj, features="percent.hb") + ggtitle("Percentage of hemoglobin genes")
    g<-g+g2
    width=width+1600
    ncol=ncol+1
  }

  if("percent.ribo" %in% colnames(obj@meta.data)){
    g3<-FeaturePlot_scCustom(obj, features="percent.ribo") + ggtitle("Percentage of ribosomal genes")
    g<-g+g3
    width=width+1600
    ncol=ncol+1
  }
  g=g+plot_layout(ncol=ncol)

  ggsave(paste0(outFile, ".genes.png"), g, width=width, height=1500, dpi=300, units="px", bg="white")
  
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
  
  cat("draw pictures ... \n")
  draw_feature_qc(prefix=paste0(outFile, ".Ident"), 
    rawobj=obj, 
    ident_name="orig.ident")

  p<-draw_dimplot(mt, paste0(outFile, ".Ident.png"), "Ident")
  if(!all(mt$Sample == mt$Ident)){
    draw_feature_qc(prefix=paste0(outFile, ".sample"), 
      rawobj=obj, 
      ident_name="sample")
    p1<-draw_dimplot(mt, paste0(outFile, ".sample.png"), "Sample")
    p<-p+p1
    width=width + nWidth * 600
  }

  if(has_batch_file){
    draw_feature_qc(prefix=paste0(outFile, ".batch"), 
      rawobj=obj, 
      ident_name="batch")
    p2<-draw_dimplot(mt, paste0(outFile, ".batch.png"), "batch")
    p<-p+p2
    width=width+nWidth * 600
  }
  
  ggsave(paste0(outFile, ".final.png"), p, width=width, height=height, dpi=300, units="px", bg="white")

  if(!is.null(qc_genes)){
    if(qc_genes != ''){
      genes<-unlist(strsplit( qc_genes, ',' ))
      g<-FeaturePlot(obj, genes, split.by="orig.ident")
      ggsave(paste0(outFile, ".qc_genes.png"), g, width=3000, height=6000, dpi=300, units="px", bg="white")
    }
  }

  if("ADT" %in% names(obj)){
    defaultAssay=DefaultAssay(obj)

    DefaultAssay(obj)="ADT"
    adt_names=rownames(obj$ADT@counts)
    writeLines(adt_names, paste0(outFile, ".ADT.txt"))
    for(adt in adt_names){
      g<-FeaturePlot(obj, features=adt, cols = c("lightgrey", "red"), min.cutoff=0.01, max.cutoff=0.99, order=TRUE) + ggtitle(paste0(adt, " protein")) + theme(aspect.ratio=1)
      ggsave(paste0(outFile, ".", adt, ".png"), g, width=5, height=5, dpi=300, units="in", bg="white")
    }

    common_genes=intersect(adt_names, rownames(obj[["RNA"]]))
    for(adt in common_genes){
      DefaultAssay(obj)="ADT"
      g1<-FeaturePlot(obj, features=adt, cols = c("lightgrey", "red"), min.cutoff=0.01, max.cutoff=0.99, order=TRUE) + ggtitle(paste0(adt, " protein")) + theme(aspect.ratio=1)
      DefaultAssay(obj)="RNA"
      g2<-FeaturePlot(obj, features=adt, cols = c("lightgrey", "red"), min.cutoff=0.01, max.cutoff=0.99, order=TRUE) + ggtitle(paste0(adt, " RNA")) + theme(aspect.ratio=1)
      g<-g1+g2+plot_layout(ncol=2)
      ggsave(paste0(outFile, ".", adt, ".common.png"), g, width=10, height=5, dpi=300, units="in", bg="white")
    }

    DefaultAssay(obj)<-defaultAssay
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

  if(!is.null(species)){
    if(tolower(species) == "mm" | tolower(species) == "mouse"){
      genes$gene = toMouseGeneSymbol(genes$gene)
    }
  }


  if(length(allgenes) > 0){
    miss_genes=setdiff(genes$gene, allgenes)
    writeLines(miss_genes, con="miss_gene.csv")

    genes<-genes[genes$gene %in% allgenes,]
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
  dd=MyGetAssayData(SCLC, assay=assay, slot="data")
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

get_dot_plot<-function(obj, group.by, gene_groups, assay="RNA", rotate.title=TRUE, use_blue_yellow_red=TRUE, dot.scale=6){
  genes=unique(unlist(gene_groups))
  assaydata=MyGetAssayData(obj, assay=assay, slot="data")
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

get_dot_width<-function(g, min_width=4000){
  group_column=if("feature.groups" %in% colnames(g$data)){"feature.groups"}else{"id"}
  if(!all(c("features.plot",group_column) %in% colnames(g$data))){
    stop(paste0("features.plot or feature.groups is not in ", paste0(colnames(g$data), collapse = ",")))
  }
  ngenes = nrow(g$data[!duplicated(g$data[,c("features.plot",group_column)]),])
  ngroups = length(unique(g$data[,group_column]))
  width=ngenes * 40 + ngroups * 30 + 400
  return(max(width, min_width))
}

get_dot_height_num<-function(ngroups, min_height=1500, height_per_entry=60, height_additional_space=1000){
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
                          species="Hs",
                          dot.scale=6){
                            
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

  genes_df=genes_df[genes_df$gene %in% rownames(obj),,drop=FALSE]
  
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

  g<-get_dot_plot(obj, group.by, gene_groups, assay, rotate.title=rotate.title, use_blue_yellow_red=use_blue_yellow_red, dot.scale=dot.scale)
  
  return(g)
}

get_sub_bubble_plot<-function(obj, obj_res, subobj, subobj_res, bubblemap_file, add_num_cell=FALSE, species=NULL, assay="RNA"){
  old_meta<-obj@meta.data
  
  obj$fake_layer=paste0("global_", unlist(obj@meta.data[,obj_res]))

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
  
  g$data<-g$data[!grepl("^global_", g$data$id),]
  if(all(!is.null(sub_levels))){
    g$data$id<-factor(g$data$id, levels=sub_levels)
  }

  return(g)
}

draw_bubble_plot<-function(obj, cur_res, cur_celltype, bubble_map_file, prefix, width=4000, height=2000, rotate.title=TRUE, species="Hs"){
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
      names(x[1:min(length(x), n_markers)])
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
  if(file_ext(obj_file) == "rds" | file_ext(obj_file) == "RDS"){
    obj=readRDS(obj_file)
  }else{
    counts=read_scrna_data(obj_file)$counts
    obj=CreateSeuratObject(counts=counts)
    if(!is.null(sample_name)){
      obj$orig.ident=sample_name
    }
  }

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
  cat("read object from ", obj_file, "\n")
  obj=read_object(obj_file, meta_rds=meta_rds, columns=columns, sample_name=sample_name)
  return(obj)
}

output_ElbowPlot<-function(obj, outFile, reduction){
  p<-ElbowPlot(obj, ndims = 40, reduction = reduction)
  ggsave(paste0(outFile, ".elbowplot.", reduction, ".png"), p, width=1500, height=1200, dpi=300, units="px", bg="white")
}

draw_feature_qc<-function(prefix, rawobj, ident_name) {
  Idents(rawobj)<-ident_name

  nsample<-length(unique(unlist(rawobj[[ident_name]])))
  
  feats<-c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo")
  if("percent.hb" %in% colnames(rawobj@meta.data)){
    feats<-c(feats, "percent.hb")
  }

  if(nsample > 50){
    ncol=1
    height=10000
  }else{
    ncol=ifelse(length(feats) > 4, 3, 2)
    height=4000
  }

  cat("draw qc voilin ...\n")
  g<-VlnPlot(rawobj, features = feats, pt.size = 0.1, ncol = ncol, raster=FALSE) + NoLegend() & theme(axis.title.x = element_blank())
  ggsave(paste0(prefix, ".qc.violin.png"), g, width=6000, height=height, dpi=300, units="px", bg="white")

  if('umap' %in% names(rawobj@reductions)){
    nfeature<-length(feats)

    by.col=nfeature>=nsample
    g<-FeaturePlot(rawobj, feats, split.by=ident_name, reduction="umap", order=T, by.col=by.col, raster=FALSE) +
      theme(aspect.ratio=1)
    if(by.col){
      width = min(50000, nsample * 700)
      height = nfeature * 700
    }else{
      width = nfeature * 700 + 300
      height = min(50000, nsample * 700)
      g = g + theme(strip.text.y=element_text(angle=0))
    }

    cat("draw qc exp ...\n")
    ggsave(paste0(prefix, ".qc.exp.png"), g, width=width, height=height, dpi=300, units="px", bg="white", limitsize = FALSE)
  }

  cat("draw qc scatter ...\n")
  p1 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE) + NoLegend()
  p2 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE) + NoLegend()
  p<-p1+p2+plot_layout(ncol=2)
  ggsave(paste0(prefix, ".qc.png"), p, width=2600, height=1200, dpi=300, units="px", bg="white")
  
  mt<-data.frame(mt=rawobj$percent.mt, Sample=unlist(rawobj[[ident_name]]), nFeature=log10(rawobj$nFeature_RNA), nCount=log10(rawobj$nCount_RNA))
  nwidth=ceiling(sqrt(nsample))
  nheight=ceiling(nsample/nwidth)

  cat("draw qc mt and read ...\n")  
  p1<-ggplot(mt, aes(y=mt,x=nCount) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ylab("Percentage of mitochondrial") + xlab("log10(number of read)") +
    facet_wrap(Sample~.) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  ggsave(paste0(prefix, ".qc.read.png"), p1, width=min(20000, 500 * nwidth + 300), height=min(10000, 500*nheight), dpi=300, units="px", bg="white", limitsize = FALSE)

  cat("draw qc mt and feature ...\n")  
  p2<-ggplot(mt, aes(y=mt,x=nFeature) ) +
    geom_bin2d(bins = 70) + 
    scale_fill_continuous(type = "viridis") + 
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    ylab("Percentage of mitochondrial") + xlab("log10(number of feature)") +
    facet_wrap(Sample~.) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
  ggsave(paste0(prefix, ".qc.feature.png"), p2, width=min(20000, 500 * nwidth + 300), height=min(10000, 500*nheight), dpi=300, units="px", bg="white", limitsize = FALSE)

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

    if(any(rawobj$project != rawobj$orig.ident)){
      meta<-rawobj@meta.data
      meta$os = paste0(meta$orig.ident, ":", meta$project)
      meta<-meta[!duplicated(meta$os),,drop=F]
      smap=split(as.character(meta$project), as.character(meta$orig.ident))
      smap2=lapply(smap,function(x){
        paste0(x, collapse=",")
      })
      ct$Project=unlist(smap2[ct$Sample])
    }

    if("batch" %in% colnames(rawobj@meta.data)){
      if(any(rawobj$batch != rawobj$orig.ident)){
        meta<-unique(rawobj@meta.data[,c("batch","orig.ident")])
        smap=split(as.character(meta$batch), as.character(meta$orig.ident))
        ct$Batch=unlist(smap[ct$Sample])
      }
    }
  }
  write.table(ct, paste0(prefix, ".cell.txt"), sep="\t", row.names=F)

  cat("draw qc cell bar...\n")  
  g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Sample)) 
  if("Batch" %in% colnames(ct)){
    g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Batch))
  }else if ("Project" %in% colnames(ct)){
    g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Project))
  }else if("Source" %in% colnames(ct)){
    g<-ggplot(ct, aes(x=Sample, y=Cell, fill=Source))
  }
  g<-g + geom_bar(stat="identity") + theme_bw3(axis.x.rotate = T) +
    theme(axis.title.x=element_blank()) 
  if(ncol(ct) == 2){
    g<-g + NoLegend()
  }
  ggsave(paste0(prefix, ".cell.bar.png"), g, width=max(3000, nrow(ct) * 60), height=2000, dpi=300, units="px", bg="white")
}

myScaleData<-function(obj, features, assay, ...){
  if(is_assay_5_plus(obj, assay)){
    if("scale.data" %in% names(obj[[assay]]@layers)){
      scaled.genes<-rownames(obj[[assay]]$scale.data)
    }else{
      scaled.genes<-c()
    }
  }else{
    scaled.genes<-rownames(obj[[assay]]@scale.data)
  }

  if(!all(features %in% scaled.genes)){
    new.genes<-unique(features, scaled.genes)
    obj=ScaleData(obj, features=new.genes, assay=assay, ... )
  }
  return(obj)
}

get_top_markers<-function(markers, n){
  markers=markers[markers$p_val_adj < 0.05,]
  top10 <- markers %>% group_by(cluster) %>% top_n(n = n, wt = .data[["avg_log2FC"]])
  return(top10)
}

get_top10_markers<-function(markers){
  return(get_top_markers(markers, 10))
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

get_dim_plot<-function(obj, group.by, label.by, label=T, title=label.by, legend.title=label.by, reduction="umap", split.by=NULL, ncol=1, random_colors=TRUE, scolors=NULL, ggplot_default_colors=FALSE, color_seed=123, ...){
  labels<-obj@meta.data[,c(group.by, label.by)]
  labels<-labels[!duplicated(labels[,group.by]),]
  labels<-labels[order(labels[,group.by]),]
  cts<-as.character(labels[,label.by])

  ngroups<-length(unlist(unique(obj[[group.by]])))
  if(all(is.null(scolors))){
    scolors = scCustomize_Palette(num_groups = ngroups, 
                                  ggplot_default_colors = ggplot_default_colors, 
                                  color_seed = color_seed)    
    #scolors = get_hue_colors(ngroups, random_colors)
  }

  g<-MyDimPlot(obj, group.by=group.by, label=label, reduction=reduction, split.by=split.by, ...)+ 
    scale_color_manual(legend.title, values=scolors, labels = cts, guide = guide_legend(ncol=ncol)) + 
    ggtitle(title)
  return(g)
}

build_dummy_cluster<-function(obj, label.by, new_cluster_name, new_cluster_name_label=paste0(new_cluster_name, "_label"), add_count=FALSE){
  if(add_count){
    label_by_count<-paste0(label.by, "_count")
    obj@meta.data=add_column_count(obj@meta.data, label.by, label_by_count)
    label.by=label_by_count
  }

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

get_dim_plot_labelby<-function(obj, label.by, title=label.by, label=T, legend.title=label.by, reduction="umap", split.by=NULL, ncol=1, label_has_cluster=FALSE, add_count=FALSE, ...){
  group.by="dummy_cluster"
  group.label="dummy_label"

  obj<-build_dummy_cluster(obj, label.by, group.by, group.label, add_count=add_count)

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

get_seurat_sum_count<-function(obj, cluster_name, min_cell_per_sample=1, target_folder="./"){
  clusterDf<-obj@meta.data
  if("seurat_clusters" %in% colnames(clusterDf)){
    cts = as.character(unique(clusterDf[order(clusterDf$seurat_clusters, decreasing = F), cluster_name]))
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
    cat("  extracted", ncol(de_obj), "cells and", nrow(de_obj), "genes\n")

    dtb<-table(de_obj$orig.ident)

    snames<-names(dtb)[dtb < min_cell_per_sample]
    if(length(snames) > 0){
      cat("those samples were excluded due to cell less than", min_cell_per_sample, ": ", paste(snames, collapse = ","), "\n")
      dtb<-dtb[dtb >= min_cell_per_sample]
      de_obj<-subset(de_obj, orig.ident %in% names(dtb))
    }

    # if("counts" %in% slotNames(de_obj@assays$RNA)){
    #   ct_count<-de_obj@assays$RNA@counts
    # }else{
    #   ct_count<-de_obj@assays$RNA@layers$counts
    #   rownames(ct_count)<-rownames(de_obj)
    #   colnames(ct_count)<-colnames(de_obj)
    # }

    ct_count<-MyGetAssayData(obj = de_obj, assay = "RNA", slot = "counts")
    groupings<-unlist(de_obj$orig.ident)
    p_count<-sumcount(ct_count, groupings)

    p_file=paste0(target_folder, prefix, ".pseudo_count.csv")
    
    write.csv(p_count, p_file)
    
    res_files<-c(res_files, file_path_as_absolute(p_file))
  }

  res_df<-data.frame("cluster"=cts, "prefix"=prefixList, "pseudo_file"=res_files)
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

  cat("draw top 20 gene figure\n")
  C <- MyGetAssayData(obj=rawobj, assay="RNA", slot="counts")
  if(ncol(C) > 100000){
    cat("Too many cells, sample 100000 cells for top20 genes visualization\n")
    C <- C[,sample(1:ncol(C), 100000)]
  }
  C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  mc <- MatrixGenerics::rowMedians(C)
  most_expressed <- order(mc, decreasing = T)[20:1]
  tm <- as.matrix(Matrix::t(C[most_expressed,]))

  png(paste0(outFile, ".top20.png"), width=3000, height=2000, res=300)
  par(mar = c(4, 8, 2, 1))
  boxplot(tm, cex = 0.1, las = 1, xlab = "% total count per cell",
          col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
  dev.off()

  has_project = ifelse("project" %in% colnames(rawobj@meta.data), any(rawobj$orig.ident != rawobj$project),FALSE)
  has_sample = ifelse("sample" %in% colnames(rawobj@meta.data), any(rawobj$orig.ident != rawobj$sample), FALSE)

  cat("draw qc by orig.ident\n")
  draw_feature_qc(outFile, rawobj, "orig.ident")

  if(has_sample){
    cat("draw qc by sample\n")
    draw_feature_qc(paste0(outFile, ".sample"), rawobj, "sample")
  }

  if(has_project){
    cat("draw qc by project\n")
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

do_PCA_Integration<-function( subobj, 
                              assay, 
                              by_sctransform, 
                              method, 
                              new.reduction, 
                              orig.reduction="pca",
                              thread=1,
                              detail_prefix=NULL,
                              ignore_variable_genes=NULL) {
  if(length(ignore_variable_genes) > 0){
    vgenes=VariableFeatures(subobj)
    VariableFeatures(subobj) <- setdiff(vgenes, ignore_variable_genes)
  }

  cat("RunPCA ... \n")
  subobj <- RunPCA(object = subobj, assay=assay, verbose=FALSE)

  if(!is.null(detail_prefix)){
    output_ElbowPlot(subobj, detail_prefix, "pca")
  }

  normalization_method = ifelse(by_sctransform, "SCT", "LogNormalize")

  if(method == "FastMNNIntegration"){
    cat("IntegrateLayers by FastMNNIntegration with thread", thread, "... \n")

    is_unix = .Platform$OS.type == "unix"
    if(is_unix){
      ncores = as.numeric(thread)
      bpparam = MulticoreParam(workers = ncores)
      register(bpparam)
    }else{
      bpparam = SerialParam()
    }

    subobj <- IntegrateLayers(
      object = subobj,
      method = FastMNNIntegration,
      orig.reduction = NULL,
      assay = assay,
      new.reduction = new.reduction,
      verbose = T,
      #for fastMNN
      batch = subobj@meta.data$batch, #batch is required for integration when only one object provided
      BPPARAM = bpparam
    )
  }else if (method == "CCAIntegration" | method == "RPCAIntegration"){
    subobj <- IntegrateLayers(
      object = subobj,
      method = method,
      orig.reduction = orig.reduction,
      assay = assay,
      new.reduction = new.reduction,
      normalization.method = normalization_method,
      verbose = T
    )
  }else{
    subobj <- IntegrateLayers(
      object = subobj,
      method = method,
      orig.reduction = orig.reduction,
      assay = assay,
      new.reduction = new.reduction,
      verbose = T
    )
  }

  #The FastMNNIntegration will create a new Seurat object with cur_array (might be SCT).
  #The default assay of new object will be set to RNA. It will cause problem when we use FindNeighbors and FindClusters.
  #So we need to set the default assay to cur_assay.
  subobj[[new.reduction]]@assay.used = assay

  return(subobj)  
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
                     reduction.name = "umap",
                     redo_fastmnn = FALSE,
                     thread=1,
                     detail_prefix=NULL,
                     ignore_variable_genes=NULL){
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
    if(curreduction == "pca"){
      #https://github.com/satijalab/seurat/issues/5244
      if (by_sctransform) {
        cat(key, "use old sctransform result\n")
        #due to very limited cell numbers in small cluster, it may cause problem to redo sctransform at individual sample level, 
        #so we will keep the old data structure
        #subobj<-do_sctransform(subobj, vars.to.regress=vars.to.regress)
      }else{
        cat(key, "redo normalization\n")
        subobj<-do_normalization( subobj, 
                                  selection.method="vst", 
                                  nfeatures=2000, 
                                  vars.to.regress=vars.to.regress, 
                                  scale.all=FALSE, 
                                  essential_genes=essential_genes,
                                  ignore_variable_genes=ignore_variable_genes)
      }

      cat(key, "RunPCA\n")
      subobj<-RunPCA(subobj, npcs=cur_npcs)
    }
    
    if(curreduction == "fastmnn"){
      if(redo_fastmnn){
        fastmnn_rds = paste0(detail_prefix, ".fastmnn.rds")
        if(file.exists(fastmnn_rds)){
          cat(key, "redo FastMNN - load cached FastMNN result ...\n")
          subobj@reductions = readRDS(fastmnn_rds)
        }else{
          cat(key, "redo FastMNN ...\n")
          if(!("batch" %in% colnames(subobj))){
            subobj$batch = subobj$orig.ident
          }

          if(assay == "RNA") {
            #Based on the following Seurat v5 integration tutorial, we need to split the RNA assay by batch first,
            #and then run the NormalizeData, FindVariableFeatures and RunPCA before integration. 
            #If we use sctranform before, we will not redo the sctranform.
            #https://satijalab.org/seurat/articles/seurat5_integration
            #When using Seurat v5 assays, we can instead keep all the data in one object, but simply split the layers. 
            cat(key, "split RNA by batch ...\n")
            subobj[["RNA"]] <- split(subobj[["RNA"]], f = subobj$batch)

            cat(key, "NormalizeData ...\n")
            subobj <- NormalizeData(subobj)

            cat(key, "FindVariableFeatures ... \n")
            subobj = FindVariableFeatures(subobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
          }

          subobj = do_PCA_Integration(subobj, 
                                      assay, 
                                      by_sctransform, 
                                      method="FastMNNIntegration", 
                                      new.reduction=curreduction, 
                                      orig.reduction="pca",
                                      thread=thread,
                                      detail_prefix=detail_prefix)

          if(assay == "RNA") {
            cat(key, "JoinLayers ... \n")
            DefaultAssay(subobj) <- "RNA"
            subobj <- JoinLayers(subobj)
          }

          cat(key, "save new fastmnn result\n")
          saveRDS(subobj@reductions, fastmnn_rds)
        }
      }
    }
  }
  cat("curreduction =", curreduction, "\n")

  cat(key, "FindNeighbors by", curreduction, "\n")
  subobj<-FindNeighbors(object=subobj, 
                        assay=assay,
                        reduction=curreduction, 
                        k.param=k_n_neighbors, 
                        dims=cur_pca_dims, 
                        verbose=FALSE)

  cat(key, "FindClusters\n")
  subobj<-FindClusters( object=subobj, 
                        random.seed=random.seed, 
                        resolution=resolutions, 
                        verbose=FALSE)

  if(do_umap){
    cat(key, "RunUMAP\n")
    cur_min_dist = 0.3
    subobj<-RunUMAP(object = subobj, 
                    min.dist = cur_min_dist, 
                    reduction=curreduction, 
                    n.neighbors=u_n_neighbors, 
                    dims=cur_pca_dims, 
                    verbose = FALSE, 
                    reduction.name=reduction.name)
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

save_umap<-function(file_prefix, obj, umap_names=names(obj@reductions$umap) ){
  umap<-FetchData(obj, umap_names)
  saveRDS(umap, paste0(file_prefix, ".rds"))
  write.csv(umap, paste0(file_prefix, ".csv"))
}

get_group_colors_from_designdata<-function(designdata){      
  designUniq<-unique(designdata[,c("Group", "DisplayGroup")])
  rownames(designUniq)<-designUniq$Group
  
  controlGroup<-designUniq["control","DisplayGroup"]
  sampleGroup<-designUniq["sample","DisplayGroup"]

  groupColors<-c("blue", "red")
  names(groupColors)<-c(controlGroup, sampleGroup)
  return(groupColors)
}

get_sig_gene_figure<-function(cell_obj, sigout, designdata, sig_gene, DE_by_cell=TRUE, is_between_cluster=FALSE, log_cpm=NULL){
  groupColors<-get_group_colors_from_designdata(designdata)
  display_group_levels<-names(groupColors)

  if(!is_between_cluster){
    ddata<-designdata[!duplicated(designdata$Sample),]

    gdismap<-unlist(split(ddata$DisplayGroup, ddata$Sample))

    cell_obj@meta.data$DisplayGroup=factor(gdismap[cell_obj$orig.ident], levels=display_group_levels)
  }

  logFC<-sigout[sig_gene, "logFC"]
  FDR<-sigout[sig_gene,"FDR"]

  stopifnot(sig_gene %in% rownames(cell_obj))

  geneexp=FetchData(cell_obj,vars=c(sig_gene))
  colnames(geneexp)<-"Gene"
  colorRange<-c(min(geneexp), max(geneexp))
  fix.sc <- scale_color_gradientn(colors=c("lightgrey", "red"), limits = colorRange)
  
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
    
    p2<-MyFeaturePlot(object = cell_obj, features=as.character(sig_gene), order=T, raster=FALSE, pt.size=0.5)
    p<-p0+p1+p2+plot_layout(design="AA
BC")
    
  }else{
    p0<-ggplot(geneexp, aes(x=Sample, y=Gene, color=Group)) + 
      geom_violin() + geom_jitter(width = 0.2) + 
      theme_bw3() + theme_rotate_x_axis_label() +
      scale_color_manual(values = groupColors) +
      xlab("") + ylab("Gene Expression")
    
    if("subumap" %in% names(cell_obj@reductions)){
      p1<-MyFeaturePlot(object = cell_obj, features=sig_gene, order=T, reduction="subumap", raster=FALSE, pt.size=0.5) + xlab("UMAP_1") + ylab("UMAP_2")
    }else{
      p1<-MyFeaturePlot(object = cell_obj, features=sig_gene, order=T, reduction="umap", raster=FALSE, pt.size=0.5)
    }
    p1$data$Group=cell_obj@meta.data[rownames(p1$data), "DisplayGroup"]
    p1<-p1+facet_grid(~Group) + theme_bw3() + ggtitle("") + theme(aspect.ratio=1)
    
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
        counts = MyGetAssayData(counts, assay="RNA", slot = "counts")
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

  if(!is_seurat_object(counts)){
    counts=ceiling(counts)
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

  g<-ggplot(md, aes(cluster, db_cell_type) ) +
    geom_tile(aes(fill = cta_score)) +
    scale_fill_gradientn(colours=colors) + theme_bw3(axis.x.rotate=TRUE) + xlab("") + ylab("")

  if(!is.na(filename)){
    ggsave(filename, g, width=width, height=height, dpi=300, units="px", bg="white")
  }else{
    print(g)
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

my_FeaturePlot_scCustom <-function(seurat_object, features, ...){
  g=FeaturePlot_scCustom( seurat_object = seurat_object, 
                          features = features,
                          colors_use = viridis_plasma_dark_high, 
                          na_color = "lightgray", ...) & theme(aspect.ratio=1)
  return(g)
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
  calc_width_per_cell=50,
  label_height=500){

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

  alltbl$Category=factor(alltbl$Category, levels=(valid_columns))

  g<-ggplot(alltbl, aes(Var2, Freq, fill=Var2)) + 
    geom_bar(width=0.5, stat = "identity") + 
    facet_grid(Var1~Category, scales = "free", space='free_x') + 
    theme_bw3(TRUE) + ylab("No. cell") + xlab("") + NoLegend() +
    theme(strip.text.y = element_text(angle = 0))

  if(!is.null(bar_file)){
    height = max(1500, length(unique(alltbl$Var1)) * calc_height_per_cluster + label_height)
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

update_rownames<-function(counts, name_func){
  df=data.frame(oldgenes=rownames(counts), newgenes=name_func(rownames(counts)))
  df=df[!duplicated(df$newgenes),,drop=FALSE]
  counts=counts[df$oldgenes,,drop=FALSE]
  rownames(counts)<-name_func(rownames(counts))
  return(counts)
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
      counts=update_rownames(counts, toMouseGeneSymbol)
    }
    if (species=="Hs") {
      counts=update_rownames(counts, toupper)
    }
    sobj = CreateSeuratObject(counts = counts, project = sample_name)
    sobj$orig.ident <- sample_name
    rm(counts)
  }
  return(sobj)
}

cbind_celltype<-function(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts){
  if(is.null(cur_layermap)){
    return(cur_cts)
  }
  layer_ids<-unlist(cur_layermap[new_cluster_ids])
  names(layer_ids) <- colnames(data_norm)
  
  oldcluster<-subobj[[cluster]][[1]]
  cur_cts$seurat_clusters=oldcluster
  cur_cts$raw_cell_type<-new_cluster_ids[oldcluster]
  cur_cts$raw_seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$raw_cell_type) 
  cur_cts$cell_type<-layer_ids[oldcluster]
  cur_cts$seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$cell_type)

  return(cur_cts)
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
                           species="Hs",
                           reduction="pca",
                           ignore_variable_genes=NULL){
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

    if(pct == "Unassigned") {
      #for iteration 1, all cells are unassigned, unless harmony has already been performed in integration, we will use "pca" for clustering
      if(reduction != "pca"){
        curreduction=reduction
      }else{
        curreduction="pca"
        if(by_harmony){
          if("harmony" %in% names(subobj@reductions)){
            curreduction="harmony"
          }
        }
      }

      cat(key, "FindNeighbors by", curreduction, "\n")
      subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)

      cat(key, "FindClusters\n")
      subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolution, verbose=FALSE)
    }else{      
      if(reduction != "pca"){
        curreduction=reduction
      }else{
        curreduction="pca"
        if(by_harmony){
          if("harmony" %in% names(subobj@reductions)){
            curreduction="harmony"
          }
        }
      }
      subobj = sub_cluster(subobj = subobj, 
                            assay =  assay, 
                            by_sctransform = by_sctransform, 
                            by_harmony = by_harmony, 
                            redo_harmony = by_harmony,
                            curreduction = curreduction, 
                            k_n_neighbors = k_n_neighbors,
                            u_n_neighbors = u_n_neighbors,
                            random.seed = random.seed,
                            resolutions = resolution,
                            cur_npcs = cur_npcs, 
                            cur_pca_dims = cur_pca_dims,
                            vars.to.regress = vars.to.regress,
                            essential_genes = essential_genes,
                            key = key,
                            detail_prefix = curprefix,
                            ignore_variable_genes = ignore_variable_genes)
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
    subobj = AddMetaData(subobj, factor_by_count(cur_cts$seurat_clusters), "seurat_clusters")
    subobj = AddMetaData(subobj, factor_by_count(cur_cts$cell_type), "cell_type")
    subobj = AddMetaData(subobj, factor_by_count(cur_cts$seurat_cell_type), "seurat_cell_type")
    subobj = AddMetaData(subobj, factor_by_count(cur_cts$raw_cell_type), "raw_cell_type")
    subobj = AddMetaData(subobj, factor_by_count(cur_cts$raw_seurat_cell_type), "raw_seurat_cell_type")
    
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
    cur_width=4000 + ceiling(length(unique(subobj$seurat_clusters)) / 20) * 1600
    ggsave(umap_celltype_file, g, width=cur_width, height=4000, units="px", dpi=300, bg="white")
    files<-rbind(files, c(previous_layer, iter_name, pct, "old_umap", umap_celltype_file))

    if(pct != "Unassigned"){
      g1<-get_dim_plot(subobj, group.by="seurat_clusters", label.by="raw_seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Raw cell type in new UMAP")
      g2<-get_dim_plot(subobj, group.by="seurat_clusters", label.by="seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat cell type in new UMAP")

      g<-g1+g2+plot_layout(nrow=1)
      umap_cluster_file = paste0(curprefix, ".", pct_str, ".new_umap.png")
      ggsave(umap_cluster_file, g, width=cur_width, height=2000, units="px", dpi=300, bg="white")
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
    ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(subobj, cluster), units="px", dpi=300, bg="white")

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
                                 species,
                                 reduction="pca",
                                 ignore_variable_genes=NULL){
  meta<-obj@meta.data
  
  previous_celltypes<-unique(meta[[previous_layer]])

  files = get_empty_files()

  if(length(previous_celltypes) > 0){
    cat("cluster and annotate cell types:", paste0(previous_celltypes, ", "))
    iter = 1
    while(TRUE){
      cat("Iteration ", iter, "\n")
      
      iter_name=paste0("iter", iter)

      previous_celltypes<-sort(previous_celltypes)
      
      curprefix = paste0(prefix, ".iter", iter)

      iter_meta_file = paste0(curprefix, ".csv")
      iter_meta_rds = paste0(curprefix, ".rds")

      if (0) { # if we want to use the cached data
        if(file.exists(iter_meta_rds)){
          cat("Found previous meta file, skip this iteration\n")
          obj@meta.data = readRDS(iter_meta_rds)
          pre_disagree<-obj@meta.data[obj@meta.data[, previous_layer] != obj@meta.data[, iter_name],,drop=F]
          if(nrow(pre_disagree) > 0){
            previous_layer = iter_name
            previous_celltypes=unique(obj@meta.data[, previous_layer])
            iter = iter + 1
            next
          }else{
            cat("No more cell type to iterate\n")
            break
          }
        }
      }

      cat("  Call iterate_celltype ...\n")

      cur_by_harmony = by_harmony
      if(iter == 1){
        if(length(previous_celltypes)==1){
          if(previous_celltypes[1] == "Unassigned"){
            cur_by_harmony=FALSE
          }
        }
      }
      
      lst<-iterate_celltype(obj, 
                            previous_celltypes, 
                            previous_layer, 
                            iter_name, 
                            cur_layermap, 
                            npcs, 
                            resolution, 
                            random.seed, 
                            by_sctransform, 
                            cur_by_harmony, 
                            curprefix, 
                            iter, 
                            vars.to.regress,
                            bubblemap_file, 
                            essential_genes,
                            species=species,
                            reduction=reduction,
                            ignore_variable_genes=ignore_variable_genes)

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
  ggsave(paste0(prefix, ".", cur_layer, ".umap.png"), g, width=3300, height=3000, units="px", dpi=300, bg="white")

  if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
    g2<-get_bubble_plot(obj, 
                        NULL, 
                        cur_layer, 
                        bubblemap_file, 
                        assay="RNA",
                        species=species)
    ggsave(paste0(prefix, ".", cur_layer, ".dot.png"), g2, width=get_dot_width(g2), height=get_dot_height(obj, cur_layer), units="px", dpi=300, bg="white")
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
                      species="Hs",
                      init_layer="layer0",
                      final_layer="layer4",
                      reduction="pca",
                      ignore_variable_genes=NULL) {
  tmp_prefix=file.path(tmp_folder, prefix)
  cur_prefix=file.path(cur_folder, prefix)

  if(0) {#for debug
    previous_layer = init_layer
    cur_layer = final_layer
    cur_layermap = layer2map
    prefix = tmp_prefix
  }

  reslist1<-layer_cluster_celltype( obj = obj,
                                    previous_layer = init_layer, 
                                    cur_layer = final_layer, 
                                    cur_layermap = layer2map, 
                                    npcs = npcs, 
                                    resolution = resolution, 
                                    random.seed = random.seed, 
                                    by_sctransform = by_sctransform, 
                                    by_harmony = by_harmony, 
                                    prefix = tmp_prefix, 
                                    vars.to.regress = vars.to.regress,
                                    bubblemap_file = bubblemap_file,
                                    essential_genes = essential_genes,
                                    species=species,
                                    reduction=reduction,
                                    ignore_variable_genes=ignore_variable_genes)
  obj=reslist1$obj
  files=reslist1$files
  rm(reslist1)

  write.csv(files, paste0(cur_prefix, ".iter_png.csv"))

  celltypes<-unique(obj@meta.data[,final_layer])
  celltypes<-celltypes[order(celltypes)]
  ctdf<-data.frame("celltype"=celltypes, "resolution"=0.01)
  write.table(ctdf, paste0(tmp_prefix, ".scDynamic.celltype_res.txt"), row.names=F, sep="\t", quote=F)

  obj<-factorize_layer(obj, final_layer)
  Idents(obj)<-final_layer

  cat("Save meta data ...\n")
  saveRDS(obj@meta.data, file.path(cur_folder, paste0(prefix, ".scDynamic.meta.rds")))
  write.csv(obj@meta.data, file.path(cur_folder, paste0(prefix, ".scDynamic.meta.csv")))

  save_umap(paste0(tmp_prefix, ".scDynamic.umap"), obj)

  if(length(celltypes) > 1){
    cat("find markers for all cell types ...\n")
    all_markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)

    cat("get top 10 markers ...\n")
    all_top10<-get_top10_markers(all_markers)
    all_top10<-unique(all_top10$gene)

    cat("scale top 10 markers ...\n")
    obj<-myScaleData(obj, features=all_top10, assay="RNA")

    cat("marker gene heatmap ...\n")
    g<-MyDoHeatMap(obj, max_cell=5000, assay="RNA", features = all_top10, group.by = final_layer, angle = 90) + NoLegend()

    width<-get_heatmap_width(length(celltypes))
    height<-get_heatmap_height(length(all_top10))
    ggsave(paste0(tmp_prefix, ".", final_layer, ".heatmap.png"), g, width=width, height=height, units="px", dpi=300, bg="white")
  }
  
  cat("output_celltype_figures...\n")
  output_celltype_figures(obj, 
    final_layer, 
    tmp_prefix, 
    bubblemap_file, 
    cell_activity_database, 
    combined_ct_source, 
    group.by="orig.ident", 
    name="sample",
    species=species)

  if(!by_individual_sample){
    #output individual sample dot plot, with global scaled average gene expression.
    obj$sample_final<-paste0(obj$orig.ident, ":", obj@meta.data[,final_layer])
    g<-get_bubble_plot( obj = obj, 
                        cur_res = NA, 
                        cur_celltype = "sample_final", 
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
      dot_file = paste0(tmp_prefix, ".", final_layer, ".", sample, ".dot.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(obj, final_layer), units="px", dpi=300, bg="white")
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
      output_celltype_figures(obj, final_layer, tmp_prefix, bubblemap_file, cell_activity_database, combined_ct_source, group.by="batch", name="batch")
    }
  }

  final_clusters=paste0(final_layer, "_clusters")
  final_raw=paste0(final_layer, "_raw")
  obj$seurat_final=paste0(obj@meta.data[,final_clusters], ": ", obj@meta.data[,final_raw])

  cur_celltype=final_layer
  final_celltypes=unique(obj@meta.data[,cur_celltype])
  for(pct in final_celltypes){
    cells=colnames(obj)[obj@meta.data[,cur_celltype] == pct]
    subobj=subset(obj, cells=cells)
    g<-get_dim_plot(subobj, group.by=final_clusters, label.by="seurat_final", reduction="umap", legend.title="")

    ggsave(paste0(tmp_prefix, ".", cur_celltype, ".", celltype_to_filename(pct), ".umap.png"), g, width=2400, height=2000, units="px", dpi=300, bg="white")
  }

  if(by_individual_sample){
    tb=data.frame("cell"=colnames(obj), "cell_type"=obj@meta.data[,final_layer], category=prefix)
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

add_column_count<-function(meta, column, target_column, order_by_count=TRUE){
  tbl=table(meta[,column])
  tbl=tbl[tbl > 0]
  if(order_by_count){
    tbl=tbl[order(tbl, decreasing=T)]
  }else{
    tbl=tbl[order(names(tbl))]
  }
  new_tbl=paste0(names(tbl), " (", tbl, ")")
  names(new_tbl)=names(tbl)
  meta[,target_column]=factor(new_tbl[as.character(meta[,column])], levels=new_tbl)
  return(meta)
}

write_gsea_rnk_by_loose_criteria<-function(dge_all, groups, design, prefix){
  #for GSEA, we need to use more genes, so we decrease the filter criteria
  keep <- filterByExpr(dge_all, group=groups, min.count=1, min.total.count=5)
  dge_loose <- dge_all[keep, , keep.lib.sizes=FALSE]
  cat("  GSEA - calcNormFactors", "\n")
  dge_loose<-calcNormFactors(dge_loose, method="TMM")

  cat("  GSEA - estimateDisp\n")
  dge_loose<-estimateDisp(dge_loose, design=design)

  cat("  GSEA - glmFit", "\n")
  fit<-glmFit(dge_loose, design=design, robust=TRUE)
  fitTest<-glmLRT(fit)
  out<-topTags(fitTest, n=Inf)
  gseaFile<-paste0(prefix, "_GSEA.rnk")
  pValuesNoZero=out$table$PValue
  rMinValue=.Machine$double.eps * .Machine$double.xmin
  if (any(pValuesNoZero < rMinValue)){
    pValuesNoZero[pValuesNoZero < rMinValue] = rMinValue
  }
  rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * (-log10(pValuesNoZero)))
  #rankout<-data.frame(gene=rownames(out), sigfvalue=sign(out$table$logFC) * (-log10(out$table$PValue)))
  rankout<-rankout[order(rankout$sigfvalue, decreasing=TRUE),]
  write.table(rankout, file=gseaFile, row.names=F, col.names=F, sep="\t", quote=F)

  return(gseaFile)
}

save_volcano_plot<-function(edgeR_out_table, 
                            prefix, 
                            useRawPvalue=0, 
                            pvalue=0.05, 
                            foldChange=2, 
                            comparisonTitle="",
                            width=7,
                            height=7,
                            extensions=c("png", "pdf")){
  library(EnhancedVolcano)  
  yname=bquote(-log[10](p~value))
  if(useRawPvalue == 1){
    pCutoffCol="PValue"
  }else{
    pCutoffCol="FDR"
  }
  p<-EnhancedVolcano(edgeR_out_table,
      lab = rownames(edgeR_out_table),
      x = 'logFC',
      y = 'PValue',
      title = comparisonTitle,
      pCutoff = pvalue,
      pCutoffCol = pCutoffCol,
      FCcutoff = log2(foldChange),
      pointSize = 3.0,
      labSize = 6.0,
      colAlpha = 1,
      subtitle = NULL) + ylab(yname) + theme(plot.title = element_text(hjust = 0.5))
  for(extension in extensions){
    volcano_file=paste0(prefix, ".volcano.", extension)
    ggsave(volcano_file, p, width=width, height=height, units="in", dpi=300, bg="white")
  }
}

is_seurat_v5_object<-function(obj){
  return(Version(obj) > '5')
}

convert_seurat_v5_to_v3<-function(obj){
  if(is_seurat_v5_object(obj)){
    obj[["RNA"]] <- as(object = obj[["RNA"]], Class = "Assay")
    obj@version=package_version("3.0.0")
  }
  return(obj)
}

get_filtered_obj<-function(subobj, filter_column, ct_meta=subobj@meta.data){
  ct_tbl = table(ct_meta[,filter_column])
  ct_tbl = ct_tbl[order(ct_tbl, decreasing=T)]
  top=names(ct_tbl)[1:min(10, length(ct_tbl))]

  ct_tbl = ct_tbl / sum(ct_tbl)
  ct_tbl = ct_tbl[ct_tbl >= 0.01]
  topperc=names(ct_tbl)

  all_cts=unique(c(top, topperc))

  cur_meta = ct_meta[as.character(ct_meta[,filter_column]) %in% all_cts,]
  cells = rownames(cur_meta)
  ct_obj = subset(subobj, cells=cells)

  ct_obj@meta.data[,filter_column] = factor_by_count(ct_obj@meta.data[,filter_column])
  return(ct_obj)
}

get_scale_color_gradient2<-function(gdata, score_name, cols) {
  if(length(cols) == 3){
    midpoint=(max(gdata[,score_name]) + min(gdata[,score_name])) /2
    result = scale_color_gradient2(low=cols[1], mid=cols[2], high=cols[3], midpoint=midpoint)
  }else{
    result = scale_color_gradient(low=cols[1], high=cols[3])
  }
  return(result)
}

get_feature_plot<-function(cur_obj, cur_feature, cols, order=FALSE) {
  g=FeaturePlot(cur_obj, 
                features = cur_feature, 
                pt.size = 0.01, 
                raster=FALSE,
                order=order) + 
    xlab("UMAP_1") + ylab("UMAP_2") + 
    theme(aspect.ratio=1,
          plot.title = element_text(hjust = 0.5))

  g=g+get_scale_color_gradient2(g$data, cur_feature, cols)
  return(g)
}

get_clustered_dotplot<-function(cur_obj, features, group.by, dot.scale=6, clustered_by_gene=TRUE, clustered_by_sample=FALSE, clustering_distance="euclidean"){
  g<-DotPlot(cur_obj, features = features, group.by=group.by,dot.scale=dot.scale) + 
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) +
    scale_colour_gradient2(low = "lightgray", high = "red")

  gdata = g$data 
  na_genes = unique(as.character(gdata$features.plot[is.na(gdata$avg.exp.scaled)]))
  gdata = gdata[!(gdata$features.plot %in% na_genes),]
  gmat = acast(data=gdata, formula = id ~ features.plot, value.var="avg.exp.scaled")

  result = list()
  if(clustered_by_gene){
    if(clustering_distance == "euclidean"){
      rc = hclust(dist(t(gmat), method = "euclidean"))
    }else{
      rc = hclust(as.dist(1 - cor(gmat, use = "pa")))
    }
    ordered_genes = c(rc$labels[rc$order], na_genes)
    g$data$features.plot = factor(as.character(g$data$features.plot), levels=ordered_genes)

    result$ordered_genes = ordered_genes
  }

  if(clustered_by_sample){
    cc = hclust(as.dist(1 - cor(t(gmat), use = "pa")))
    ordered_samples = cc$labels[cc$order]
    g$data$id = factor(g$data$id, levels=ordered_samples)

    result$ordered_samples = ordered_samples
  }

  result$g = g
  return(result)
}

get_freq_table<-function(df, column){
  df |> 
    dplyr::count(!!sym(column)) |> 
    dplyr::arrange(desc(n)) |>
    dplyr::rename("count"="n")
}

get_category_with_min_percentage<-function(obj, category, min_perc){
  tbl=get_freq_table(obj@meta.data, category)
  min_cells=sum(tbl$count) * min_perc
  tbl=tbl[tbl$count>min_cells,]
  major_cells=colnames(obj)[obj@meta.data[,category] %in% tbl[,category]]
  result_obj<-subset(obj, cells=major_cells)
  return(result_obj)
}

process_move=function(move_formula, move_parts, cur_meta){
  if(length(move_parts) != 4){
    stop(paste0("    move formula should be MOVE:cluster:column:celltypes:to_cluster, now we get: ", move_formula))
  }

  move_cluster=move_parts[1]
  cur_name=unique(cur_meta$cur_layer[cur_meta$seurat_clusters_str==move_cluster])
  move_column=move_parts[2]
  move_celltypes_str=move_parts[3]
  move_celltypes=unlist(strsplit(move_celltypes_str, ","))
  to_cluster=suppressWarnings(as.numeric(move_parts[4]))

  if(is.na(to_cluster)){
    to_seurat_clusters = FALSE
    to_cluster=move_parts[4]
  }else{
    to_seurat_clusters = TRUE
  }

  cat("    moving", move_celltypes_str, "annotated by", move_column, "from cluster", move_cluster, "to", to_cluster, "\n")
  if(!move_column %in% colnames(cur_meta)){
    stop(paste0("column ", move_column, " not exists"))
  }

  if(move_cluster == "-1" | move_cluster == "*"){
    is_move = cur_meta[,move_column] %in% move_celltypes
  }else{
    is_move = cur_meta$seurat_clusters_str==move_cluster & cur_meta[,move_column] %in% move_celltypes
  }

  # even the cells are moved before, we still need to move them again if it is not deleted
  #final_move = is_move & (!cur_meta$is_moved)
  final_move = is_move & cur_meta$seurat_clusters[is_move] >= -1

  move_cells=sum(final_move)

  if(to_seurat_clusters){
    cur_meta$seurat_clusters[final_move] = to_cluster
    cur_meta$seurat_clusters_str<-as.character(cur_meta$seurat_clusters)        
  }else{
    cur_meta$cur_layer[final_move] = to_cluster
    cur_meta$is_moved[final_move] = TRUE
  }
  cat("       cells moved:", move_cells, "\n")

  return(cur_meta)
}

process_filter=function(filter_action, filter_formula, filter_parts, cur_meta){
  if(length(filter_parts) != 3){
    stop(paste0("  filter formula should be DELETE_or_KEEP:cluster:column:celltypes, now we get: ", filter_formula))
  }
  
  filter_cluster=filter_parts[1]
  filter_column=filter_parts[2]
  filter_celltypes_str=filter_parts[3]
  filter_celltypes=unlist(strsplit(filter_celltypes_str, ","))

  if(!filter_action %in% c("DELETE", "KEEP")){
    stop(paste0("  filter action should be DELETE or KEEP, now we get: ", filter_action, "for", filter_formula))
  }

  if(filter_column == "*" & filter_action == "DELETE"){
    #delete all cells in the cluster
    delete_clusters=strsplit(filter_cluster, ",")[[1]]
    cur_name=unique(cur_meta$cur_layer[cur_meta$seurat_clusters_str %in% delete_clusters])
    cat("    deleting all cells in cluster", filter_cluster, "(", paste0(cur_name, collapse=","), ")\n")
    is_deleted = cur_meta$seurat_clusters_str %in% delete_clusters
  }else{
    if (filter_cluster == -1){
      cat("   ", filter_action, filter_celltypes_str, "annotated by", filter_column, "from all sub-clusters\n")
      if(filter_action == "KEEP"){
        is_deleted = !(cur_meta[,filter_column] %in% filter_celltypes)
      }else{
        is_deleted = cur_meta[,filter_column] %in% filter_celltypes
      }
    }else{
      cur_name=unique(cur_meta$cur_layer[cur_meta$seurat_clusters_str==filter_cluster])
      cat("   ", filter_action, filter_celltypes_str, "annotated by", filter_column, "from cluster", filter_cluster, "(", cur_name, ")\n")
      if(!filter_column %in% colnames(cur_meta)){
        stop(paste0("column ", filter_column, " not exists"))
      }

      if(filter_action == "KEEP"){
        is_deleted = cur_meta$seurat_clusters_str==filter_cluster & !(cur_meta[,filter_column] %in% filter_celltypes)
      }else{
        is_deleted = cur_meta$seurat_clusters_str==filter_cluster & cur_meta[,filter_column] %in% filter_celltypes
      }
    }
  }

  #if the cells are moved, we don't delete them. If the cells are deleted, we don't count them again
  final_deleted = is_deleted & (!cur_meta$is_moved) & cur_meta$seurat_clusters != -10000

  delete_cells=sum(final_deleted)
  cur_meta$seurat_clusters[final_deleted] = -10000
  cur_meta$seurat_clusters_str<-as.character(cur_meta$seurat_clusters)
  cat("      cells deleted:", delete_cells, "\n")

  return(cur_meta)
}

process_merge=function(action_formula, action_parts, cur_meta) {
  if(length(action_parts) != 2){
    stop(paste0("  rename formula should be MERGE:from_clusters:to_cluster, now we get: ", action_formula))
  }
  
  from_clusters_str=action_parts[1]
  from_clusters=unlist(strsplit(from_clusters_str, ","))
  to_cluster=as.numeric(action_parts[2])
  if(is.na(to_cluster)){
    stop(paste0("wrong to_cluster value ", action_parts[2]))
  }

  cat("    merge cluster:", from_clusters_str, "to", to_cluster, "\n")

  is_merge = cur_meta$seurat_clusters_str %in% from_clusters

  # if the cells are moved, we don't merge them
  final_merge = is_merge & (!cur_meta$is_moved)

  cur_meta[final_merge, "seurat_clusters"] = to_cluster
  cur_meta$seurat_clusters_str<-as.character(cur_meta$seurat_clusters)
  cat("      cells merged:", sum(final_merge), "\n")

  return(cur_meta)
}

process_rename=function(action_formula, 
                        action_parts, 
                        cur_meta, 
                        condition_column) {
  if(length(action_parts) != 2){
    stop(paste0("  rename formula should be RENAME:cluster:name, now we get", action_formula))
  }
  
  rename_cluster=action_parts[1]
  rename_clusters=unlist(strsplit(rename_cluster, ","))
  rename_name=action_parts[2]

  if(all(rename_clusters == "-1" | rename_clusters == "*")){
    is_rename = !cur_meta$is_moved
  } else {
    is_rename = (cur_meta[,condition_column] %in% rename_clusters) & !cur_meta$is_moved
  }

  rename_cells=sum(is_rename)

  cur_meta[is_rename, "cur_layer"] = rename_name

  cat("   cells renamed:", rename_cells, "\n")

  return(cur_meta)
}

process_actions=function(ct_tbl, cur_meta, condition_column="seurat_clusters_str"){
  if(nrow(cur_meta) == 0){
    stop("cur_meta is empty")
  }
  cur_meta$seurat_clusters_str<-as.character(cur_meta$seurat_clusters)
  if(nrow(ct_tbl) > 0){
    cur_meta$is_moved=FALSE
    action_formula= ct_tbl$V1[1]
    for(action_formula in ct_tbl$V1){
      cat("  action formula:", action_formula, "\n")
      action_parts=unlist(strsplit(action_formula, ":"))
      action_type=action_parts[1]
      action_parts=action_parts[2:length(action_parts)]

      if(action_type == "MOVE"){
        cur_meta=process_move(move_formula=action_formula, 
                              move_parts=action_parts, 
                              cur_meta=cur_meta)
      }else if(action_type == "KEEP" | action_type == "DELETE"){
        cur_meta=process_filter(filter_action=action_type, 
                                filter_formula=action_formula, 
                                filter_parts=action_parts, 
                                cur_meta=cur_meta)
      }else if(action_type == "MERGE"){
        cur_meta=process_merge( action_formula, 
                                action_parts, 
                                cur_meta)
      }else if(action_type == "RENAME"){
        cur_meta=process_rename(action_formula, 
                                action_parts, 
                                cur_meta, 
                                condition_column)
      }else{
        stop(paste0("wrong action type:", action_type))
      }
    }
  }

  cur_meta$cur_layer[cur_meta$seurat_clusters < -1] = "DELETE"
  return(cur_meta)
}

get_nhood_stat_ct <- function(comp_milo, annotation_column){
  nhoods_sce = nhoods(comp_milo)
  # assign cell types for nhoods 
  nhood_stat_ct = data.frame(Nhood = 1:ncol(nhoods_sce) , Nhood_center = colnames(nhoods_sce))
  nhood_stat_ct = miloR::annotateNhoods(comp_milo , 
                                        nhood_stat_ct , 
                                        coldata_col = annotation_column)
  return(nhood_stat_ct)
}
