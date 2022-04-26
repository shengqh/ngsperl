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

min_markers<-20

previous_layer<-myoptions$celltype_layer
cur_layer<-myoptions$output_layer
seurat_cur_layer=paste0("seurat_", cur_layer)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

pca_dims<-1:npcs

tiers<-read.table(myoptions$HLA_panglao5_file, sep="\t", header=T)

remove_subtype_of=remove_subtype
cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file, curated_markers_file=myoptions$curated_markers_file)

prefix<-outFile

if(!exists("obj")){
  obj=readRDS(parFile1)
  if(is.list(obj)){
    obj<-obj$obj
  }
  obj@meta.data = readRDS(parFile2)
  Idents(obj)<-previous_layer
}

if(parSampleFile2 != ""){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  obj<-obj[!(rownames(obj) %in% ignore_genes),]
}

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  bubble_genes<-unique(genes_df$`Marker Gene`)
}

resolutions=seq(from = 0.1, to = 0.9, by = 0.1)

layer4map<-split(tiers$Layer4, tiers$Celltype.name)

meta = obj@meta.data
curprefix = prefix

assay=ifelse(by_sctransform, "SCT", "RNA")
clusters=paste0(assay, "_snn_res.", resolutions)

celltypes<-unlist(meta[[previous_layer]])
tblct<-table(celltypes)
tblct<-tblct[order(tblct, decreasing = T)]

previous_celltypes<-names(tblct)

allmarkers<-NULL
allcts<-NULL
cluster_index=0
pct<-previous_celltypes[5]
for(pct in previous_celltypes){
  key = paste0(previous_layer, ": ", pct, ":")
  cells<-rownames(meta)[meta[,previous_layer] == pct]
  subobj<-subset(obj, cells=cells)
  
  stopifnot(all(subobj[[previous_layer]] == pct))
  
  pca_npcs<-min(round(length(cells)/2), 50)
  
  cur_npcs=min(pca_npcs, npcs)
  cur_pca_dims=1:cur_npcs
  
  curprefix<-paste0(prefix, ".", previous_layer, "_", gsub(" ", "_", pct))
  
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
  
  cat(key, "Find best resolution\n")
  cluster = clusters[1]
  lastCluster = clusters[1]
  for(cluster in clusters){
    cat("  ", cluster, "\n")
    Idents(subobj)<-cluster
    if(length(unique(Idents(subobj))) == 1){
      lastCluster = cluster
      next
    }
    
    markers=FindAllMarkers(subobj, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    markers=markers[markers$p_val_adj < 0.05,]
    nmarkers=unlist(lapply(unique(Idents(subobj)), function(x){sum(markers$cluster==x)}))
    if(all(nmarkers >= min_markers)){
      lastCluster = cluster
    }else{
      break
    }
  }
  
  Idents(subobj)<-lastCluster
  has_one_cluster<-length(unique(Idents(subobj))) == 1
  if(has_one_cluster){
    cts<-subobj[[previous_layer]]
    cts$resolution<-0.1
    cts$seurat_clusters<-cluster_index
    cts[,cur_layer] = cts[,previous_layer]
    cts[,seurat_cur_layer] = paste0(cts$seurat_clusters, ": ", cts[,cur_layer])
    cluster_index=cluster_index+1
    markers<-FindMarkers(obj, ident.1 = pct, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    markers=markers[markers$p_val_adj < 0.05,]
    top10 <- markers %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
    top10genes<-unique(top10$gene)
    allmarkers<-c(allmarkers, top10genes)
  }else{
    markers=FindAllMarkers(subobj, only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    markers=markers[markers$p_val_adj < 0.05,]
    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
    top10genes<-unique(top10$gene)
    
    allmarkers<-c(allmarkers, top10genes)
    
    data.norm=get_seurat_average_expression(subobj, lastCluster)
    predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
    layer_ids<-names(predict_celltype$max_cta)
    names(layer_ids) <- colnames(data.norm)
    
    oldcluster<-subobj[[lastCluster]][[1]]
    newct<-layer_ids[oldcluster]
    newdf<-data.frame("layer"=newct)
    colnames(newdf)<-cur_layer

    cts<-subobj[[previous_layer]]
    cur_resolution=gsub(".+_snn_res.", "", lastCluster)
    cts$resolution<-cur_resolution
    
    subclusters<-sort(unique(subobj@meta.data[[lastCluster]]))
    newclusters<-cluster_index+as.numeric(as.character(subclusters))
    names(newclusters)<-subclusters
    cts$seurat_clusters<-unlist(newclusters[oldcluster])
    cts<-cbind(cts, newdf)
    cts[,seurat_cur_layer] = paste0(cts$seurat_clusters, ": ", cts[,cur_layer])

    ct<-cts[,c("seurat_clusters", seurat_cur_layer)]
    ct<-unique(ct)
    ct<-ct[order(ct$seurat_clusters),]
    
    subobj<-AddMetaData(subobj, cts$seurat_clusters, col.name = "seurat_clusters")
    subobj<-AddMetaData(subobj, cts[,cur_layer], col.name = cur_layer)
    subobj<-AddMetaData(subobj, factor(cts[,seurat_cur_layer], levels=ct[,seurat_cur_layer]), col.name = seurat_cur_layer)
    
    cluster_index<-max(cts$seurat_clusters)+1
    
    cat("RUNUMAP\n")
    subobj<-RunUMAP(object = subobj, reduction=curreduction, dims=cur_pca_dims, verbose = FALSE)
    
    g<-DimPlot(subobj, group.by = "seurat_clusters", label=T) + ggtitle(paste0(pct, ": res", cur_resolution) ) + scale_color_discrete(labels = ct[,seurat_cur_layer])
    if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
      layout <- "ABB"
      g2<-get_bubble_plot(subobj, "seurat_clusters", cur_layer, bubblemap_file, assay)
      g<-g+g2+plot_layout(design=layout)
      width=6300
    }else{
      width=2300
    }
    
    png(paste0(curprefix, ".umap.png"), width=width, height=2000, res=300)
    print(g)
    dev.off()
    
    width<-max(3000, min(10000, length(unique(subobj$seurat_clusters)) * 150 + 1000))
    height<-max(3000, min(10000, length(top10genes) * 60 + 1000))
    
    png(paste0(curprefix, ".top10.heatmap.png"), width=3000, height=3000, res=300)
    g<-DoHeatmap(subobj, features = top10genes, group.by = seurat_cur_layer, angle = 90) + NoLegend()
    print(g)
    dev.off()
  }
  
  allcts<-rbind(allcts, cts)
}

write.csv(allcts, paste0(outFile, ".", previous_layer, "_to_", cur_layer, ".cts.csv"))

stopifnot(nrow(allcts) == ncol(obj))
allcts<-allcts[colnames(obj),]

ct<-allcts[,c("seurat_clusters", seurat_cur_layer)]
ct<-unique(ct)
ct<-ct[order(ct$seurat_clusters),]

obj<-AddMetaData(obj, allcts$seurat_clusters, col.name = "seurat_clusters")
obj<-AddMetaData(obj, allcts[,cur_layer], col.name = cur_layer)
obj<-AddMetaData(obj, factor(allcts[,seurat_cur_layer], levels=ct[,seurat_cur_layer]), col.name = seurat_cur_layer)

allmarkers<-unique(allmarkers)

width<-max(3000, min(10000, length(unique(obj$seurat_clusters)) * 150 + 1000))
height<-max(3000, min(20000, length(allmarkers) * 60 + 1000))

png(paste0(prefix, ".top10.heatmap.png"), width=width, height=height, res=300)
g<-DoHeatmap(obj, features = allmarkers, group.by = seurat_cur_layer, angle = 90) + NoLegend()
print(g)
dev.off()

g<-DimPlot(obj, group.by = "seurat_clusters", label=T) + ggtitle(cur_layer)+
      scale_color_discrete(labels = ct[,seurat_cur_layer])
if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
  g<-g+get_bubble_plot(obj, "seurat_clusters", cur_layer, bubblemap_file, assay)
  g<-g+plot_layout(ncol = 2, widths = c(4, 6))
  width=11000
}else{
  width=4300
}

png(paste0(prefix, ".umap.png"), width=width, height=4000, res=300)
print(g)
dev.off()

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))
