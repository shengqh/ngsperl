rm(list=ls()) 
outFile='P9061'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_merge/result/P9061.final.rds'
parFile2=''
parFile3='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/essential_genes/result/P9061.txt'


setwd('/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_merge_dynamic_01_call/result')

### Parameter setting end ###

source("scRNA_func.r")
library(plyr)
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

by_sctransform<-is_one(myoptions$by_sctransform)
reduction<-myoptions$reduction
npcs<-as.numeric(myoptions$pca_dims)

species=myoptions$species
markerfile<-myoptions$db_markers_file
remove_subtype<-myoptions$remove_subtype
annotate_tcell<-is_one(myoptions$annotate_tcell)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(by_sctransform, "SCT", "RNA")
by_harmony<-reduction=="harmony"
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)
redo_harmony<-is_one(myoptions$redo_harmony, 0)
resolution=as.numeric(myoptions$dynamic_by_one_resolution)

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
cell_activity_database$cellType=cell_activity_database$cellType[names(cell_activity_database$cellType) %in% tiers$Celltype.name]
tiers<-tiers[tiers$Celltype.name %in% names(cell_activity_database$cellType),]

cts<-read.table(parSampleFile4, header=F, sep="\t", stringsAsFactors = F)
combined_ct<-unlist(split(cts$V2, cts$V1))
layers = paste0("Layer", c(1:4))
for(layer in layers){
  tiers[,layer] = mapvalues(tiers[,layer], names(combined_ct), combined_ct)
}

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

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes, species = myoptions$species)
  bubble_genes<-unique(genes_df$gene)
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

cbind_celltype<-function(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts){
  if(is.null(cur_layermap)){
    return(cur_cts)
  }
  layer_ids<-unlist(cur_layermap[new_cluster_ids])
  names(layer_ids) <- colnames(data_norm)
  
  oldcluster<-subobj[[cluster]][[1]]
  cur_cts$seurat_clusters=oldcluster
  cur_cts$cell_type<-layer_ids[oldcluster]
  cur_cts$seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$cell_type)

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

if(0){
  previous_layer<-"layer2"
  cur_layer="layer3"
  previous_layermap=layer2map
  cur_layermap=layer3map
  previous_celltypes<-unique(obj@meta.data[[previous_layer]])
  iter=1
  check_pre_layer=TRUE
}

if(1){
  previous_layer<-"layer0"
  cur_layer="layer1"
  previous_layermap=NA
  cur_layermap=layer1map
  previous_celltypes<-unique(obj@meta.data[[previous_layer]])
  iter=1
  check_pre_layer=FALSE
}

iterate_celltype<-function(obj, previous_celltypes, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, iter, check_pre_layer, vars.to.regress){
  meta = obj@meta.data
  curprefix = paste0(prefix, ".", previous_layer, ".iter", iter)
  
  assay=ifelse(by_sctransform, "SCT", "RNA")
  
  all_pre_cts<-NULL
  all_cur_cts<-NULL
  pct<-previous_celltypes[length(previous_celltypes)]
  
  #previous_celltypes<-c("Platelets")
  for(pct in previous_celltypes){
    pct_str = gsub('[/\ ]', "_", pct)

    key = paste0(previous_layer, ": iter", iter, ": ", pct, ":")
    cells<-rownames(meta)[meta[,previous_layer] == pct]
    if(length(cells) == 0){#no cell left for this cell type
      next
    }
    
    #g0<-DimPlot(obj, label=F, cells.highlight =cells) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))
    g0<-DimPlot(obj, label=F, cells.highlight =cells) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))

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
      subobj = sub_cluster(subobj, 
                            assay, 
                            by_sctransform, 
                            by_harmony, 
                            redo_harmony, 
                            curreduction, 
                            k_n_neighbors,
                            u_n_neighbors,
                            random.seed,
                            resolution,
                            cur_npcs, 
                            cur_pca_dims,
                            vars.to.regress, 
                            essential_genes, 
                            key,
                            do_umap = TRUE,
                            reduction.name = "umap")
    }else{
      cat(key, "FindNeighbors\n")
      subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)

      cat(key, "FindClusters\n")
      subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolution, verbose=FALSE)
    }
    
    cat(key, "Cell type annotation\n")
    pre_cts<-subobj[[previous_layer]]
    cur_cts<-subobj[[previous_layer]]

    cluster<-"seurat_clusters"
    data_norm=get_seurat_average_expression(subobj, cluster)
    
    predict_celltype<-ORA_celltype(data_norm,cell_activity_database$cellType,cell_activity_database$weight)
    
    new_cluster_ids<-names(predict_celltype$max_cta)
    
    cur_cts<-cbind_celltype(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts)
    subobj = AddMetaData(subobj, cur_cts$cell_type, "layer")
    
    if(check_pre_layer){
      pre_cts<-cbind_celltype(subobj, data_norm, cluster, new_cluster_ids, previous_layermap, pre_cts)
      subobj = AddMetaData(subobj, pre_cts$cell_type, "pre_layer")
    }
    
    #using RNA assay for visualization
    DefaultAssay(subobj)<-assay

    if(check_pre_layer){
      g1<-DimPlot(subobj, reduction="oumap", group.by = "pre_layer", label=T) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle(paste0(pct, ": ", previous_layer, ": post"))
    }else{
      g1<-DimPlot(subobj, group.by = previous_layer, label=T) + ggtitle(paste0(pct, ": ", previous_layer, ": pre"))
    }
    g2<-DimPlot(subobj, group.by = "layer", label=T) + ggtitle(paste0(pct, ": ", cur_layer))

    if(check_pre_layer){
      g3<-DimPlot(subobj, group.by = "pre_layer", label=T) + ggtitle(paste0(pct, ": ", previous_layer, ": post"))
      g<-g0+g1+g3+g2
      width=4500
      height=4000
      ncol=2
    }else{
      g<-g0+g1+g2
      width=6750
      height=2000
      ncol=3
    }
    png(paste0(curprefix, ".", pct_str, ".umap.png"), width=width, height=height, res=300)
    print(g)
    dev.off()

    dot_width=4400
    g<-get_bubble_plot(subobj, cur_res=cluster, "layer", bubblemap_file, assay="RNA", orderby_cluster = FALSE)
    png(paste0(curprefix, ".", pct_str, ".dot.png"), width=dot_width, height=get_dot_height(subobj, cluster), res=300)
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

layer_cluster_celltype<-function(obj, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, check_pre_layer, vars.to.regress, unique_celltype){
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
      
      lst<-iterate_celltype(obj, previous_celltypes, previous_layer, previous_layermap, cur_layer, cur_layermap, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, iter, check_pre_layer, vars.to.regress)
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
    g2<-get_bubble_plot(obj, NA, cur_layer, bubblemap_file, assay)
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

obj<-layer_cluster_celltype(obj, "layer0", NULL, "layer1", layer1map, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, FALSE, vars.to.regress, NULL)
#obj@meta.data=readRDS(paste0(outFile, ".layer0_to_layer1.meta.rds"))

obj<-layer_cluster_celltype(obj, "layer1", layer1map, "layer2", layer2map, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, TRUE, vars.to.regress, NULL)
#obj@meta.data=readRDS(paste0(outFile,".layer1_to_layer2.meta.rds"))

obj<-layer_cluster_celltype(obj, "layer2", layer2map, "layer3", layer3map, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, TRUE, vars.to.regress, NULL)
#obj@meta.data=readRDS(paste0(outFile,".layer2_to_layer3.meta.rds"))

obj<-layer_cluster_celltype(obj, "layer3", layer3map, "layer4", layer4map, npcs, resolution, random.seed, by_sctransform, by_harmony, prefix, TRUE, vars.to.regress, unique_celltype)
#obj@meta.data=readRDS(paste0(outFile,".layer3_to_layer4.meta.rds"))

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

obj<-myScaleData(obj, all_top10, "RNA")
if(ncol(obj) > 5000){
  subsampled <- obj[, sample(colnames(obj), size=5000, replace=F)]
  g<-DoHeatmap(subsampled, assay="RNA", group.by="layer4", features=all_top10)
  rm(subsampled)
}else{
  g<-DoHeatmap(obj, assay="RNA", group.by="layer4", features=all_top10)
}

width<-max(3000, min(10000, length(unique(Idents(obj))) * 150 + 1000))
height<-max(3000, min(10000, length(all_top10) * 60 + 1000))
png(paste0(outFile, ".layer4.heatmap.png"), width=width, height=height, res=300)
print(g)
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

g<-get_celltype_marker_bubble_plot( obj = obj, 
                                    group.by = "layer4", 
                                    cellType = cell_activity_database$cellType,
                                    weight = cell_activity_database$weight,
                                    n_markers = 5, 
                                    combined_ct=combined_ct)

png(paste0(outFile, ".layer4.ct_markers.bubbleplot.png"), width=5500, height=3000, res=300)
print(g)
dev.off()

if("batch" %in% colnames(obj@meta.data)){
  output_barplot(obj, "batch", "layer4", paste0(outFile, ".batch_cluster.png"))
}

save_highlight_cell_plot(paste0(prefix, ".", "layer4", ".cell.png"), obj, group.by = "layer4", reduction = "umap");

library('rmarkdown')
rmarkdown::render("seurat_scDynamic_one_resolution.rmd",output_file=paste0(outFile,".html"))

