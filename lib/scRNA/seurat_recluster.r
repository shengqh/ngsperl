
source("scRNA_func.r")

library(Seurat)
library(ggplot2)
library(patchwork)
library(kableExtra)
library(dplyr)

random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_Mt_rRNA= ifelse(myoptions$Remove_Mt_rRNA == "FALSE", FALSE, TRUE)
resolution=as.numeric(myoptions$resolution)
species=myoptions$species
markerfile<-myoptions$markers_file

by_integration<-ifelse(myoptions$by_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

prefix<-myoptions$prefix

pca_dims<-1:as.numeric(myoptions$pca_dims)

cts_cluster=read_cell_cluster_file(parFile2)

recluster_celltypes<-myoptions$recluster_celltypes
if(recluster_celltypes != ""){
  recluster_celltypes<-unlist(strsplit(recluster_celltypes, ";"))
  cts=cts_cluster[,myoptions$celltype_name]
  all_celltypes=unique(cts_cluster[,myoptions$celltype_name])
  missed=recluster_celltypes[!(recluster_celltypes %in% all_celltypes)]
  if(length(missed) > 0){
    stop(paste0("There are missed celltypes ", paste0(missed, collapse="/"), " not in file ", parFile2))
  }
}

celltype_cells=data.frame(Cell=rownames(cts_cluster), CellType=cts_cluster[,myoptions$celltype_name])
cclist=split(celltype_cells$Cell, celltype_cells$CellType)

finalList<-readRDS(parFile1)
obj<-finalList$obj

obj[[myoptions$cluster_name]]=cts_cluster[colnames(obj),myoptions$cluster_name]

seurat_colors<-finalList$seurat_colors
seurat_cellactivity_colors<-finalList$seurat_cellactivity_colors

allrds<-NULL

ct<-recluster_celltypes[1]
for(ct in recluster_celltypes){
  ctPrefix<-paste0(prefix, ".", gsub(' ', '_', ct))
  cells=unlist(cclist[ct])
  cluster_obj<-subset(obj, cells=cells)

  dir.create(ctPrefix, showWarnings = FALSE)
  ctPrefix=file.path(ctPrefix, ctPrefix)
  
  png(paste0(ctPrefix, ".pre.png"), width=3000, height=3000, res=300)
  p<-DimPlot(object = cluster_obj, reduction = 'umap', label=TRUE, group.by=myoptions$cluster_name) + ggtitle("")
  print(p)
  dev.off()
  
  rdsFile = paste0(ctPrefix, ".rds")
  if(!file.exists(rdsFile)){
    if(by_sctransform){
      cluster_obj <- SCTransform(cluster_obj, verbose = FALSE)
    }

    obj_markers <- run_cluster(cluster_obj, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, resolution, random.seed)
    cluster_obj<-obj_markers$object

    clusterMarkers<-obj_markers$markers %>% group_by(cluster)
    write.csv(clusterMarkers, file=paste0(ctPrefix, ".allmarkers.csv"), row.names=F, quote = F)

    if ("Seurat" %in% names(sessionInfo()$otherPkgs) & grepl("^4",sessionInfo()$otherPkgs$Seurat$Version)) { #Seurat version4
      logFcColName="avg_log2FC"
    } else {
      logFcColName="avg_logFC"
    }

    top10 <- obj_markers$markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[[logFcColName]])
    top10marker_file = paste0(ctPrefix, ".top10markers.csv")
    write.csv(top10, file=top10marker_file, row.names=F, quote = F)

    cur_display_markers=rownames(top10)
    dot_filename=paste0(ctPrefix, ".top10markers.dot.pdf")
    pdf(file=dot_filename, width=14, height=7)
    g=DotPlot(obj, features=cur_display_markers, assay="SCT", group.by="seurat_cellactivity_clusters" ) + 
      xlab("") + ylab("") + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, hjust=1))
    print(g)
    dev.off()
    
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

  allrds <-rbind(allrds, data.frame(category=ct, rds=rdsFile))

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

write.csv(allrds, file=paste0(outFile, ".rds.list"), row.names=F)
