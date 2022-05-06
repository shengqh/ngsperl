library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(digest)
library(heatmap3)
library(cowplot)
library(scales)
library(stringr)
library(htmltools)
library(patchwork)
library(glmGamPoi)

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
reduction<-myoptions$reduction
npcs<-as.numeric(myoptions$pca_dims)

species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")
by_harmony<-reduction=="harmony"
regress_by_percent_mt<-ifelse(myoptions$regress_by_percent_mt == "1", TRUE, FALSE)

min_markers=as.numeric(myoptions$best_resolution_min_markers)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

tiers<-read.table(myoptions$HLA_panglao5_file, sep="\t", header=T)

remove_subtype_of=remove_subtype
cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file, curated_markers_file=myoptions$curated_markers_file)

prefix<-outFile

obj<-readRDS(parFile1)
if(is.list(obj)){
  obj<-obj$obj
}

if(parSampleFile2 != ""){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  cat("removing genes in", ignore_gene_files$V1, "\n")
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  obj<-obj[!(rownames(obj) %in% ignore_genes),]
}

obj<-subset(obj, cells=colnames(obj)[1:1000])

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  bubble_genes<-unique(genes_df$`Marker Gene`)
}

resolutions=seq(from = 0.1, to = 0.9, by = 0.1)

layer1map<-split(tiers$Layer1, tiers$Celltype.name)
layer2map<-split(tiers$Layer2, tiers$Celltype.name)
layer3map<-split(tiers$Layer3, tiers$Celltype.name)
layer4map<-split(tiers$Layer4, tiers$Celltype.name)

obj[["layer0"]]<-"Unassigned"

cbind_celltype<-function(subobj, data.norm, cluster, new.cluster.ids, cur_layermap, cur_cts){
  if(is.null(cur_layermap)){
    return(cur_cts)
  }
  layer_ids<-unlist(cur_layermap[new.cluster.ids])
  names(layer_ids) <- colnames(data.norm)
  
  oldcluster<-subobj[[cluster]][[1]]
  newct<-layer_ids[oldcluster]
  newdf<-data.frame("cell_type"=newct, "seurat_clusters"=oldcluster)
  newdf$seurat_cell_type<-paste0(newdf$ct, ": ", newdf$cl)
  cur_cts<-cbind(cur_cts, newdf)
  return(cur_cts)
}

previous_layer<-"layer0"
cur_layer="layer1"
previous_layermap=NULL
cur_layermap=layer1map
previous_celltypes<-unique(obj@meta.data[[previous_layer]])
iter=1

iterate_celltype<-function(obj, previous_celltypes, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, iter, check_pre_layer, vars.to.regress){
  meta = obj@meta.data
  curprefix = paste0(prefix, ".", previous_layer, ".iter", iter)
  
  assay=ifelse(by_sctransform, "SCT", "RNA")
  clusters=paste0(assay, "_snn_res.", resolutions)
  
  all_pre_cts<-NULL
  all_cur_cts<-NULL
  pct<-previous_celltypes[length(previous_celltypes)]
  
  #previous_celltypes<-c("Platelets")
  for(pct in previous_celltypes){
    key = paste0(previous_layer, ": iter", iter, ": ", pct, ":")
    cells<-rownames(meta)[meta[,previous_layer] == pct]
    subobj<-subset(obj, cells=cells)
    
    stopifnot(all(subobj[[previous_layer]] == pct))
    
    pca_npcs<-min(round(length(cells)/2), 50)
    
    cur_npcs=min(pca_npcs, npcs)
    cur_pca_dims=1:cur_npcs
    

    if(by_harmony){
      cat(key, "harmony\n")
      DefaultAssay(subobj)<-"RNA"
      subobj[['harmony']]<-NULL
      subobj<-do_harmony(subobj, cur_npcs, batch_file="")
      curreduction="harmony"
    }else if (by_sctransform) {
      cat(key, "sctransform\n")
      DefaultAssay(subobj)<-"RNA"
      subobj[['SCT']]<-NULL
      subobj<-SCTransform(subobj, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes = F, verbose = FALSE)
      subobj<-RunPCA(subobj, npcs=pca_npcs)
      curreduction="pca"
    }else{
      cat(key, "normalization\n")
      DefaultAssay(subobj)<-"RNA"
      subobj<-NormalizeData(subobj)
      subobj<-FindVariableFeatures(subobj)

      var.genes<-VariableFeatures(subobj)
      var.genes<-unique(c(var.genes, essential_genes))

      subobj<-ScaleData(subobj, vars.to.regress = vars.to.regress, features = var.genes)
      subobj<-RunPCA(subobj, npcs=pca_npcs)
      curreduction="pca"
    }
    cat(key, "FindClusters\n")
    subobj<-FindNeighbors(object=subobj, reduction=curreduction, dims=cur_pca_dims, verbose=FALSE)
    subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolutions, verbose=FALSE)
    
    cat(key, "RunUMAP\n")
    subobj<-RunUMAP(object = subobj, reduction=curreduction, dims=cur_pca_dims, verbose = FALSE)
    
    cat(key, "Find best resolution\n")
    pre_cts<-subobj[[previous_layer]]
    cur_cts<-subobj[[previous_layer]]

    cluster<-find_best_resolution(subobj, clusters, min.pct, logfc.threshold, min_markers)
    Idents(subobj)<-cluster

    cat(key, "Cell annotation for best resolution", cluster, "\n")
    data.norm=get_seurat_average_expression(subobj, cluster)
    predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
    new.cluster.ids<-names(predict_celltype$max_cta)
    cur_cts<-cbind_celltype(subobj, data.norm, cluster, new.cluster.ids, cur_layermap, cur_cts)
      
    if(check_pre_layer){
      pre_cts<-cbind_celltype(subobj, data.norm, cluster, new.cluster.ids, previous_layermap, pre_cts)
    }
    
    subobj = AddMetaData(subobj, cur_cts$cell_type, "layer")

    if(check_pre_layer){
      subobj = AddMetaData(subobj, pre_cts$cell_type, "pre_layer")
    }
    
    g1<-DimPlot(subobj, group.by = previous_layer, label=T) + ggtitle(paste0(pct, ": ", previous_layer, ": pre"))
    g2<-DimPlot(subobj, group.by = "layer", label=T) + ggtitle(paste0(pct, ": ", cur_layer))

    if(check_pre_layer){
      g3<-DimPlot(subobj, group.by = "pre_layer", label=T) + ggtitle(paste0(pct, ": ", previous_layer, ": post"))
      g<-g1+g3+g2
      width=6900
    }else{
      g<-g1+g2
      width=4600
    }

    #saveRDS(subobj, paste0(curprefix, "_", gsub(" ", "_", pct), ".rds"))

    if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
      if(check_pre_layer){
        layout <- "
ABC
DDD
"
      }else{
        layout <- "
AB
CC
"
      }
      g4<-get_bubble_plot(subobj, NA, "layer", bubblemap_file, assay="RNA")
      g<-g+g4+plot_layout(design=layout)
      height=4000
    }else{
      if(check_pre_layer){
        g<-g+plot_layout(ncol=3)
      }else{
        g<-g+plot_layout(ncol=2)
      }
      height=2000
    }
    png(paste0(curprefix, "_", gsub(" ", "_", pct), ".png"), width=width, height=height, res=300)
    print(g)
    dev.off()

    all_pre_cts<-rbind(all_pre_cts, pre_cts)
    all_cur_cts<-rbind(all_cur_cts, cur_cts)
    
    if(previous_layer == "layer0"){
      obj[['umap']] = subobj[['umap']]
    }
    rm(subobj)
  }
  return(list("all_pre_cts"=all_pre_cts, "all_cur_cts"=all_cur_cts))
}

layer_cluster_celltype<-function(obj, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, check_pre_layer, vars.to.regress){
  meta<-obj@meta.data
  
  previous_celltypes<-unique(meta[[previous_layer]])
  
  iter = 1
  while(TRUE){
    cat("Iteration ", iter, "\n")
    
    iter_meta_file = paste0(prefix, ".", previous_layer, ".iter", iter, ".csv")

    lst<-iterate_celltype(obj, previous_celltypes, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, iter, check_pre_layer, vars.to.regress)
    all_pre_cts<-lst$all_pre_cts
    all_cur_cts<-lst$all_cur_cts
    
    stopifnot(all(rownames(all_cur_cts) %in% colnames(obj)))
    obj@meta.data[rownames(all_cur_cts), cur_layer]=unlist(all_cur_cts[, "cell_type"])
    
    if(check_pre_layer){
      stopifnot(all(rownames(all_pre_cts) %in% colnames(obj)))
      pre_disagree<-all_pre_cts[all_pre_cts[, previous_layer] != all_pre_cts[,'cell_type'],,drop=F]
      if(nrow(pre_disagree) > 0){
        cat("Found unmatched cell type\n")
        print(table(pre_disagree[,previous_layer], pre_disagree[,'cell_type']))
        
        obj@meta.data[rownames(pre_disagree), previous_layer] = unlist(pre_disagree$layer)
        write.csv(obj@meta.data, iter_meta_file)
        previous_celltypes = unique(c(unlist(pre_disagree[,previous_layer]), unlist(pre_disagree$layer)))
      }else{
        write.csv(obj@meta.data, iter_meta_file)
        break
      }
      iter = iter + 1
    }else{
      write.csv(obj@meta.data, iter_meta_file)
      break
    }
  }
  
  g<-DimPlot(obj, group.by = cur_layer, label=T)
  if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
    g2<-get_bubble_plot(obj, NA, cur_layer, bubblemap_file, assay="RNA")
    g<-g+g2+plot_layout(ncol = 2, widths = c(3, 5))
    width=6400
  }else{
    width=2400
  }
  png(paste0(prefix, ".", previous_layer, "_to_", cur_layer, ".png"), width=width, height=2000, res=300)
  print(g)
  dev.off()

  write.csv(obj@meta.data, file=paste0(prefix, ".", previous_layer, "_to_", cur_layer, ".meta.csv"))

  write.csv(all_cur_cts, file=paste0(prefix, ".", previous_layer, "_to_", cur_layer, ".details.csv"))
  
  return(obj)
}

obj<-layer_cluster_celltype(obj, "layer0", NULL, "layer1", layer1map, npcs, resolutions, random.seed, by_sctransform, FALSE, prefix, FALSE, vars.to.regress)

obj<-layer_cluster_celltype(obj, "layer1", layer1map, "layer2", layer2map, npcs, resolutions, random.seed, by_sctransform, FALSE, prefix, TRUE, vars.to.regress)

obj<-layer_cluster_celltype(obj, "layer2", layer2map, "layer3", layer3map, npcs, resolutions, random.seed, by_sctransform, FALSE, prefix, TRUE, vars.to.regress)

obj<-layer_cluster_celltype(obj, "layer3", layer3map, "layer4", layer4map, npcs, resolutions, random.seed, by_sctransform, FALSE, prefix, TRUE, vars.to.regress)

saveRDS(obj@meta.data, paste0(prefix, ".scDynamic.meta.rds"))
write.csv(obj@meta.data, paste0(prefix, ".scDynamic.meta.csv"))
