
#source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(patchwork)

run_cluster<-function(object, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, resolution, random.seed){
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
  object <- RunUMAP(object = object, dims=pca_dims, verbose = FALSE)
  object <- FindNeighbors(object = object, dims=pca_dims, verbose=FALSE)
  object <- FindClusters(object=object, verbose=FALSE, random.seed=random.seed, resolution=resolution)
  
  if (by_sctransform) {
    markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }else{
    markers <- FindAllMarkers(object, features=var.genes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }
  markers <- markers[markers$p_val_adj < 0.01,]
  return(list(object=object, markers=markers))
}

ORA_celltype<-function(medianexp,cellType,weight){
  ORA_result<-matrix(NA, nrow=length(cellType),ncol=dim(medianexp)[2])
  CTA_result<-matrix(0,nrow=length(cellType),ncol=dim(medianexp)[2])
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
      CTA_result[i,j]<-sum(exp_ss*weight_ss)/(length(exp_ss)^0.3)
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

random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
options<-split(options_table$V1, options_table$V2)

Mtpattern= options$Mtpattern
rRNApattern=options$rRNApattern
Remove_Mt_rRNA= ifelse(options$Remove_Mt_rRNA == "FALSE", FALSE, TRUE)
resolution=as.numeric(options$resolution)
species=options$species
markerfile<-options$markers_file

by_integration<-ifelse(options$by_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(options$by_sctransform == "0", FALSE, TRUE)

prefix<-options$prefix

pca_dims<-1:as.numeric(options$pca_dims)

recluster_celltypes<-options$recluster_celltypes
recluster_celltypes<-unlist(strsplit(recluster_celltypes, ";"))

cts_cluster=read_cell_cluster_file(parFile2)
cts=cts_cluster[,options$celltype_name]
all_celltypes=unique(cts_cluster[,options$celltype_name])
missed=recluster_celltypes[!(recluster_celltypes %in% all_celltypes)]
if(length(missed) > 0){
  stop(paste0("There are missed celltypes ", paste0(missed, collapse="/"), " not in file ", parFile2))
}

celltype_cells=data.frame(Cell=rownames(cts_cluster), CellType=cts_cluster[,options$celltype_name])
cclist=split(celltype_cells$Cell, celltype_cells$CellType)

finalList<-readRDS(parFile1)
obj<-finalList$obj

obj[[options$cluster_name]]=cts_cluster[colnames(obj),options$cluster_name]

seurat_colors<-finalList$seurat_colors
seurat_cellactivity_colors<-finalList$seurat_cellactivity_colors

ct<-recluster_celltypes[1]
for(ct in recluster_celltypes){
  ctPrefix<-paste0(prefix, ".", gsub(' ', '_', ct))
  cells=unlist(cclist[ct])
  cluster_obj<-subset(obj, cells=cells)
  
  png(paste0(ctPrefix, ".pre.png"), width=3000, height=3000, res=300)
  p<-DimPlot(object = cluster_obj, reduction = 'umap', label=TRUE, group.by=options$cluster_name) + ggtitle("")
  print(p)
  dev.off()
  
  rdsFile = paste0(ctPrefix, ".rds")
  if(!file.exists(rdsFile)){
    cluster_obj <- SCTransform(cluster_obj, verbose = FALSE)
    obj_markers <- run_cluster(cluster_obj, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, resolution, random.seed)
    cluster_obj<-obj_markers$object
    
    clusters<-cluster_obj@active.ident
    sumcounts<-t(apply(GetAssayData(cluster_obj,assay="RNA",slot="counts"),1,function(x){tapply(x,clusters,sum)}))
    logsumcounts<-log2(sumcounts+1)
    data.quantileAll <- apply(logsumcounts, 2, function(x){quantile(x, 0.75)})
    
    norm_method=""
    if(any(data.quantileAll == 0)){
      norm_method = ".normByTotal"
      data.all <- apply(logsumcounts, 2, sum)
      data.all<-data.all / median(data.all)
      data.norm <- t(t(logsumcounts) / data.all)
    }else{
      norm_method = ".normByUpQuantile"
      data.quantileAll<-data.quantileAll / median(data.quantileAll)
      data.norm <- t(t(logsumcounts) / data.quantileAll)
    }
    
    colnames(sumcounts)<-paste0("Cluster", colnames(sumcounts))
    write.csv(sumcounts, file=paste0(ctPrefix, ".cluster.count.csv"))
    
    oldname<-colnames(data.norm)
    colnames(data.norm)<-paste0("Cluster", oldname)
    write.csv(data.norm, file=paste0(ctPrefix, ".cluster", norm_method, ".csv"))
    colnames(data.norm)<-oldname
    
    predict_celltype<-ORA_celltype(data.norm,finalList$cell_activity_database$cellType,finalList$cell_activity_database$weight)
    
    new.cluster.ids<-names(predict_celltype$max_cta)
    seurat_clusters<-unlist(cluster_obj[["seurat_clusters"]])
    names(new.cluster.ids) <- levels(seurat_clusters)
    cluster_obj[["cellactivity_clusters"]] <- new.cluster.ids[unlist(cluster_obj[["seurat_clusters"]])]
    
    clusterDf<-data.frame(seurat=unlist(cluster_obj[["seurat_clusters"]]), cellactivity=unlist(cluster_obj[["cellactivity_clusters"]]))
    clusterDf$seurat_cellactivity<-paste0(clusterDf$seurat, " : ", clusterDf$cellactivity)
    seurat_cellactivity<-clusterDf$seurat_cellactivity

    clusterDf<-clusterDf[order(clusterDf$seurat),]
    seurat_cellactivity<-factor(seurat_cellactivity, levels=unique(clusterDf$seurat_cellactivity))
    cluster_obj[["seurat_cellactivity_clusters"]] <- seurat_cellactivity
    
    saveRDS(cluster_obj, file=paste0(ctPrefix, ".rds"))
  }else{
    cluster_obj<-readRDS(rdsFile)
  }
  png(paste0(ctPrefix, ".post.png"), width=3000, height=3000, res=300)
  p<-DimPlot(object = cluster_obj, reduction = 'umap', label=TRUE, group.by="seurat_cellactivity_clusters") + ggtitle("")
  print(p)
  dev.off()
  
  samples<-unique(unlist(cluster_obj[["orig.ident"]]))
  nc=ceiling(sqrt(length(samples)))
  nn=ceiling(length(samples) / nc)
  
  png(paste0(ctPrefix, ".post_sample.png"), width= nc*1000+200, height=nn*1000, res=300)
  p<-DimPlot(object = cluster_obj, reduction = 'umap', label=F, group.by="seurat_cellactivity_clusters", split.by="orig.ident", combine=F)
  p<-p[[1]] + facet_wrap(~orig.ident) + ggtitle("")
  print(p)
  dev.off()
}
