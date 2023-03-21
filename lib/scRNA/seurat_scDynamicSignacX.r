rm(list=ls()) 
outFile='AG3669'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1='C:/projects/data/h_gelbard_lab/projects/20220609_scRNA_3669_sct/seurat_sct_merge/result/AG3669.final.rds'
parFile2=''
parFile3='C:/projects/data/h_gelbard_lab/projects/20220609_scRNA_3669_sct/essential_genes/result/AG3669.txt'
parFile4='C:/projects/data/h_gelbard_lab/projects/20220609_scRNA_3669_sct/seurat_sct_merge_SignacX/result/AG3669.meta.rds'


setwd('C:/projects/data/h_gelbard_lab/projects/20220609_scRNA_3669_sct/seurat_sct_merge_01_dynamic_signacX/result')

### Parameter setting end ###

source("scRNA_func.r")
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
remove_subtype_of<-myoptions$remove_subtype
annotate_tcell<-ifelse(myoptions$annotate_tcell == "0", FALSE, TRUE)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")
by_harmony<-reduction=="harmony"
regress_by_percent_mt<-ifelse(myoptions$regress_by_percent_mt == "1", TRUE, FALSE)

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

cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file, curated_markers_file=myoptions$curated_markers_file)
cell_activity_database$cellType=cell_activity_database$cellType[names(cell_activity_database$cellType) %in% tiers$Celltype.name]
tiers<-tiers[tiers$Celltype.name %in% names(cell_activity_database$cellType),]
rownames(tiers)<-tiers$Celltype.name

cmap<-data.frame("signacx"=c("B.memory", "B.naive",
                             "DC", "Mon.Classical", "Mon.NonClassical", "Neutrophils", "Monocytes", "Macrophages", 
                             "NK", "T.CD4.naive", "T.CD4.memory", "T.regs", "T.CD8.naive", "T.CD8.memory", "T.CD8.cm","T.CD8.em", 
                             "Endothelial", "Fibroblasts", "Epithelial", "Unclassified"),
                 "panglao"=c("B cells", "B cells", 
                             "Dendritic cells", "Monocytes", "Monocytes", "Neutrophils","Monocytes", "Macrophages", 
                             "NK cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells",
                             "Endothelial cells", "Fibroblasts", "Epithelial cells", "Unclassified"))

smap<-unlist(split(cmap$panglao, cmap$signacx))

#write.csv(tiers, "tiers.csv")

prefix<-outFile

if(!exists('obj')){
  obj<-readRDS(parFile1)
  if(is.list(obj)){
    obj<-obj$obj
  }
}

if(parSampleFile2 != ""){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  cat("removing genes in", ignore_gene_files$V1, "\n")
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  obj<-obj[!(rownames(obj) %in% ignore_genes),]
}

if(parSampleFile3 != ""){
  umap_min_dist = read.table(parSampleFile3, sep="\t")
  umap_min_dist_map = unlist(split(umap_min_dist$V1, umap_min_dist$V2))
}else{
  umap_min_dist_map = c("layer0" = 0.3,"layer1" = 0.3,"layer2" = 0.3,"layer3" = 0.3,"layer4" = 0.3)
}

signacx<-readRDS(parFile4)

#write.csv(cmap, "cmap.csv")
if(!all(signacx$signacx_CellStates %in% names(smap))){
  cs<-signacx$signacx_CellStates[!(signacx$signacx_CellStates %in% names(smap))]
  stop("not all signacx in map: ", paste0(unique(cs), collapse = ", "))
}

signacx$signacx <- unlist(lapply(signacx$signacx_CellStates, function(x){
  smap[as.character(x)]
}))

obj$signacx<-signacx$signacx

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes, species = myoptions$species)
  bubble_genes<-unique(genes_df$gene)
}

if(myoptions$dynamic_by_one_resolution != ""){
  resolutions=c(as.numeric(myoptions$dynamic_by_one_resolution))
}else{
  resolutions=seq(from = 0.1, to = 0.9, by = 0.1)
}

#find cell type equals in layer2/layer3/layer4, for example, the T cells after removing sub celltypes
get_unique234_celltype<-function(tiers){
  layer_name = "Layer2"
  map<-split(tiers$Celltype.name, tiers[,layer_name])
  unique_map = Filter(function(x)(length(x) == 1), map)
  #unique_map2 = unique_map[Filter(function(x)(names(unique_map)[x] == unique_map[x][1]), c(1:length(unique_map)))]
  result=unlist(unique_map)
  names(result) = names(unique_map)
  return(result)
}

unique_celltype<-get_unique234_celltype(tiers)

layer1map<-split(tiers$Layer1, tiers$Celltype.name)
layer2map<-split(tiers$Layer2, tiers$Celltype.name)
layer3map<-split(tiers$Layer3, tiers$Celltype.name)
layer4map<-split(tiers$Layer4, tiers$Celltype.name)

obj[["layer0"]]<-"Unassigned"

cbind_celltype<-function(subobj, cluster, new.cluster.ids, cur_layermap, cur_cts){
  if(is.null(cur_layermap)){
    return(cur_cts)
  }
  layer_ids<-unlist(cur_layermap[new.cluster.ids])
  names(layer_ids) <- names(new.cluster.ids)
  
  oldcluster<-subobj[[cluster]][[1]]
  newct<-layer_ids[oldcluster]
  newdf<-data.frame("ct"=newct, "cl"=oldcluster)
  newdf$ct_cl<-paste0(newdf$ct, ": ", newdf$cl)
  colnames(newdf)<-paste0(cluster, "_", colnames(newdf))
  cur_cts<-cbind(cur_cts, newdf)
  return(cur_cts)
}

vote_celltype<-function(cts, assay){
  cls<-cts[,grepl(paste0(assay, "_snn_res.+_ct$"), colnames(cts)), drop=F]
  
  votes<-apply(cls, 1, function(x){
    tb<-table(unlist(x))
    tbwhich<-which(tb==max(tb))
    return(names(tb)[tbwhich[1]])
  })
  cts$layer<-unlist(votes)
  
  return(cts)
}

previous_layer<-"layer0"
cur_layer="layer1"
previous_layermap=NA
cur_layermap=layer1map
previous_celltypes<-unique(obj@meta.data[[previous_layer]])
check_pre_layer=FALSE
iter=1

iterate_celltype<-function(obj, previous_celltypes, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, iter, check_pre_layer, vars.to.regress, umap_min_dist_map){
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
    if(length(cells) == 0){#no cell left for this cell type
      next
    }
    
    subobj<-subset(obj, cells=cells)
    
    stopifnot(all(subobj[[previous_layer]] == pct))
    
    pca_npcs<-min(round(length(cells)/2), 50)
    cur_npcs=min(pca_npcs, npcs)
    cur_pca_dims=1:cur_npcs

    k_n_neighbors<-min(cur_npcs, 20)
    u_n_neighbors<-min(cur_npcs, 30)

    DefaultAssay(subobj)<-assay
    
    #for dynamic clustering, we will not redo normalization/sctranform
    if(by_harmony){
      cat(key, "harmony\n")
      #https://github.com/satijalab/seurat/issues/5289
      #So we can either (1) re-do the complete workflow after subsetting including the re-integration with harmony or 
      #(2) we can only re-run FindNeighbors and FindClusters. Is this correct? In our case, (2) worked better.
      curreduction="harmony"
    }else if (by_sctransform) {
      #https://github.com/satijalab/seurat/issues/2812
      #One thing to add. When you subset the clustering, Re-run PCA actually has more impacts than re-run SCTransform
      cat(key, "sctransform/PCA\n")
      #subobj[['SCT']]<-NULL
      #subobj<-SCTransform(subobj, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes = F, verbose = FALSE)
      subobj<-RunPCA(subobj, npcs=cur_npcs)
      curreduction="pca"
    }else{
      cat(key, "normalization/PCA\n")
      # subobj<-NormalizeData(subobj)
      # subobj<-FindVariableFeatures(subobj)

      # var.genes<-VariableFeatures(subobj)
      # var.genes<-unique(c(var.genes, essential_genes))

      # subobj<-ScaleData(subobj, vars.to.regress = vars.to.regress, features = var.genes)
      subobj<-RunPCA(subobj, npcs=pca_npcs)
      curreduction="pca"
    }
    
    cat(key, "FindClusters\n")
    subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)
    subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolutions, verbose=FALSE)
    
    if(previous_layer != "layer0") {
      cat(key, "RunUMAP\n")
      cur_min_dist = umap_min_dist_map[previous_layer]
      subobj<-RunUMAP(object = subobj, min.dist = cur_min_dist, reduction=curreduction, n.neighbors=u_n_neighbors, dims=cur_pca_dims, verbose = FALSE)
    }
    
    cat(key, "Cell type annotation\n")
    pre_cts<-subobj[[previous_layer]]
    cur_cts<-subobj[[previous_layer]]
    cluster = clusters[1]
    for(cluster in clusters){
      cat("  ", cluster, "\n")
      
      predict_ct<-table(unlist(subobj[[cluster]]), subobj$signacx)
      
      new.cluster.ids<-colnames(predict_ct)[apply(predict_ct,1,which.max)]
      names(new.cluster.ids) = rownames(predict_ct)
      
      is_unclassified <- which(new.cluster.ids == "Unclassified")
      if(length(is_unclassified) > 0){
        data.norm=get_seurat_average_expression(subobj, cluster)
        predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
        ora_cts<-names(predict_celltype$max_cta)
        new.cluster.ids[is_unclassified] = ora_cts[is_unclassified] 
      }
      
      
      cur_cts<-cbind_celltype(subobj, cluster, new.cluster.ids, cur_layermap, cur_cts)
      
      if(check_pre_layer){
        pre_cts<-cbind_celltype(subobj, cluster, new.cluster.ids, previous_layermap, pre_cts)
      }
    }
    
    cur_cts<-vote_celltype(cur_cts, assay)
    stopifnot(nrow(cur_cts) == ncol(subobj))
    stopifnot(all(rownames(cur_cts) == colnames(subobj)))

    subobj = AddMetaData(subobj, cur_cts$layer, "layer")

    if(check_pre_layer){
      pre_cts<-vote_celltype(pre_cts, assay)
      stopifnot(nrow(pre_cts) == ncol(subobj))
      stopifnot(all(rownames(pre_cts) == colnames(subobj)))

      subobj = AddMetaData(subobj, pre_cts$layer, "pre_layer")
    }
    
    #using RNA assay for visualization
    DefaultAssay(subobj)<-"RNA"

    cls<-cur_cts[,grepl(paste0(assay, "_snn_res.+_ct$"), colnames(cur_cts)), drop=F]
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

    g<-NULL
    for (cl in colnames(cls)){
      ct<-paste0(cl, "_cl")
      subobj[[ct]]=unlist(cur_cts[,ct])

      meta_cl<-gsub("_ct$", "", cl)

      cur_meta<-subobj@meta.data
      cur_meta<-cur_meta[!duplicated(cur_meta[,meta_cl]),]
      cur_meta<-cur_meta[order(cur_meta[,meta_cl]),]

      gg<-DimPlot(subobj, group.by = meta_cl, label=T) + ggtitle(cl) +
        scale_color_discrete(labels = cur_meta[,ct])
      if(is.null(g)){
        g=gg
      }else{
        g<-g+gg
      }
    }
    if(length(resolutions) == 1){
      png(paste0(curprefix, "_", gsub(" ", "_", pct), ".detail.png"), width=3300, height=3000, res=300)
      print(g)
      dev.off()
    }else{
      png(paste0(curprefix, "_", gsub(" ", "_", pct), ".detail.png"), width=7000, height=4800, res=300)
      g<-g+plot_layout(ncol=3)
      print(g)
      dev.off()
    }

    stopifnot(!any(rownames(cur_cts) %in% rownames(all_cur_cts)))

    all_pre_cts<-rbind(all_pre_cts, pre_cts)
    all_cur_cts<-rbind(all_cur_cts, cur_cts)
    
    if(previous_layer == "layer0"){
      obj[['umap']] = subobj[['umap']]
    }
    rm(subobj)
  }
  return(list("all_pre_cts"=all_pre_cts, "all_cur_cts"=all_cur_cts))
}

layer_cluster_celltype<-function(obj, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, check_pre_layer, vars.to.regress, unique_celltype, umap_min_dist_map){
  meta<-obj@meta.data
  
  previous_celltypes<-unique(meta[[previous_layer]])
  
  if(!is.null(unique_celltype)){
    ucs<-previous_celltypes[previous_celltypes %in% names(unique_celltype)]
    for(uc in ucs){
      obj@meta.data[meta[, previous_layer] == uc, cur_layer]=unique_celltype[uc]
    }
    previous_celltypes=previous_celltypes[!(previous_celltypes %in% ucs)]
  }
  
  if(length(previous_celltypes) > 0){
    cat("cluster and annotate cell types:", paste0(previous_celltypes, ", "))
    iter = 1
    while(TRUE){
      cat("Iteration ", iter, "\n")
      previous_celltypes<-previous_celltypes[order(previous_celltypes)]
      
      iter_meta_file = paste0(prefix, ".", previous_layer, ".iter", iter, ".csv")
      iter_meta_rds = paste0(prefix, ".", previous_layer, ".iter", iter, ".rds")
      
      lst<-iterate_celltype(obj, previous_celltypes, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, iter, check_pre_layer, vars.to.regress, umap_min_dist_map)
      all_pre_cts<-lst$all_pre_cts
      all_cur_cts<-lst$all_cur_cts
      
      stopifnot(all(rownames(all_cur_cts) %in% colnames(obj)))
      obj@meta.data[rownames(all_cur_cts), cur_layer]=unlist(all_cur_cts[, "layer"])
      
      if(check_pre_layer){
        stopifnot(all(rownames(all_pre_cts) %in% colnames(obj)))
        pre_disagree<-all_pre_cts[all_pre_cts[, previous_layer] != all_pre_cts[,'layer'],,drop=F]
        if(nrow(pre_disagree) > 0){
          cat("Found unmatched cell type\n")
          print(table(pre_disagree[,previous_layer], pre_disagree[,'layer']))
          
          obj@meta.data[rownames(pre_disagree), previous_layer] = unlist(pre_disagree$layer)
          write.csv(obj@meta.data, iter_meta_file)
          saveRDS(obj@meta.data, iter_meta_rds)
          previous_celltypes = unique(c(unlist(pre_disagree[,previous_layer]), unlist(pre_disagree$layer)))
        }else{
          write.csv(obj@meta.data, iter_meta_file)
          saveRDS(obj@meta.data, iter_meta_rds)
          break
        }
        iter = iter + 1
      }else{
        write.csv(obj@meta.data, iter_meta_file)
        saveRDS(obj@meta.data, iter_meta_rds)
        break
      }
    }
  }
  
  #using RNA assay for visualization
  DefaultAssay(obj)<-"RNA"

  g<-DimPlot(obj, group.by = cur_layer, label=T)

  png(paste0(prefix, ".", cur_layer, ".final.png"), width=3300, height=3000, res=300)
  print(g)
  dev.off()

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
  saveRDS(obj@meta.data, file=paste0(prefix, ".", previous_layer, "_to_", cur_layer, ".meta.rds"))

  write.csv(all_cur_cts, file=paste0(prefix, ".", previous_layer, "_to_", cur_layer, ".details.csv"))
  
  return(obj)
}

obj<-layer_cluster_celltype(obj, "layer0", NULL, "layer1", layer1map, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, FALSE, vars.to.regress, NULL, umap_min_dist_map)

obj<-layer_cluster_celltype(obj, "layer1", layer1map, "layer2", layer2map, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, TRUE, vars.to.regress, NULL, umap_min_dist_map)

obj<-layer_cluster_celltype(obj, "layer2", layer2map, "layer3", layer3map, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, TRUE, vars.to.regress, NULL, umap_min_dist_map)
#obj@meta.data=readRDS("iSGS_airway.layer2_to_layer3.meta.rds")
obj<-layer_cluster_celltype(obj, "layer3", layer3map, "layer4", layer4map, npcs, resolutions, random.seed, by_sctransform, by_harmony, prefix, TRUE, vars.to.regress, unique_celltype, umap_min_dist_map)

celltypes<-unique(obj$layer4)
celltypes<-celltypes[order(celltypes)]
ctdf<-data.frame("celltype"=celltypes, "resolution"=0.01)
write.table(ctdf, paste0(prefix, ".scDynamic.celltype_res.txt"), row.names=F, sep="\t", quote=F)

obj<-factorize_layer(obj, "layer4")
Idents(obj)<-"layer4"

saveRDS(obj@meta.data, paste0(prefix, ".scDynamic.meta.rds"))
write.csv(obj@meta.data, paste0(prefix, ".scDynamic.meta.csv"))

#find markers for all cell types
all_markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
all_top10<-get_top10_markers(all_markers)
all_top10<-unique(all_top10$gene)

width<-max(3000, min(10000, length(unique(Idents(obj))) * 150 + 1000))
height<-max(3000, min(10000, length(all_top10) * 60 + 1000))

obj<-myScaleData(obj, all_top10, "RNA")
png(paste0(outFile, ".layer4.heatmap.png"), width=width, height=height, res=300)
DoHeatmap(obj, assay="RNA", group.by="layer4", features=all_top10)
dev.off()

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

output_barplot(obj, "orig.ident", "layer4", paste0(outFile, ".ident_cluster.png"))

if("batch" %in% colnames(obj@meta.data)){
  output_barplot(obj, "batch", "layer4", paste0(outFile, ".batch_cluster.png"))
}
