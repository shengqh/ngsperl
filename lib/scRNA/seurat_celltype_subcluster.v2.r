library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(digest)
library(cowplot)
library(scales)
library(stringr)
library(htmltools)
library(patchwork)
library(data.table)

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
assay=ifelse(by_sctransform, "SCT", "RNA")
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
  obj<-read_object(parFile1, parFile2)
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

resolutions=c(seq(from = 0.01, to = 0.09, by = 0.02), seq(from = 0.1, to = 0.5, by = 0.1))

layer4map<-split(tiers$Layer4, tiers$Celltype.name)

meta = obj@meta.data
curprefix = prefix

clusters_prefix = paste0(assay, "_snn_res.", format(resolutions, nsmall = 2))
clusters=paste0(assay, "_snn_res.", resolutions)
names(clusters_prefix) = clusters

celltypes<-unlist(meta[[previous_layer]])
tblct<-table(celltypes)
tblct<-tblct[order(tblct, decreasing = T)]

previous_celltypes<-names(tblct)
#previous_celltypes<-c("B cells")

DefaultAssay(obj)<-assay

allmarkers<-NULL
allcts<-NULL
cluster_index=0
pct<-previous_celltypes[1]
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
    #due to very limited cell numbers in small cluster, it may cause problem to redo sctransform and harmony, 
    #so we will keep the old data structure
    #subobj<-do_harmony(subobj, by_sctransform, regress_by_percent_mt, FALSE, "", pca_dims)
    curreduction="harmony"
  }else{
    if (by_sctransform) {
      cat(key, "sctransform\n")
      #due to very limited cell numbers in small cluster, it may cause problem to redo sctransform at individual sample level, 
      #so we will keep the old data structure
      #subobj<-do_sctransform(subobj, vars.to.regress=vars.to.regress)
    }else{
      cat(key, "normalization\n")
      subobj<-do_normalization(subobj, selection.method="vst", nfeatures=3000, vars.to.regress=vars.to.regress, scale.all=FALSE, essential_genes=essential_genes)
    }
    subobj<-RunPCA(subobj, npcs=pca_npcs)
    curreduction="pca"
  }

  cat(key, "FindClusters\n")
  subobj<-FindNeighbors(object=subobj, reduction=curreduction, dims=cur_pca_dims, verbose=FALSE)
  subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolutions, verbose=FALSE)

  cat(key, "RUNUMAP\n")
  subobj<-RunUMAP(object = subobj, assay=assay, reduction=curreduction, dims=cur_pca_dims, verbose = FALSE)

  cat(key, "Find marker genes\n")
  cluster = clusters[1]
  markers_map = list()
  for(cluster in clusters){
    cat("  ", cluster, "\n")

    Idents(subobj)<-cluster
    if(length(unique(Idents(subobj))) == 1){#only 1 cluster
      cat("    only one cluster, pass\n")
      next
    }
    
    cluster_prefix = paste0(curprefix, ".", clusters_prefix[cluster])
    cur_resolution=gsub(".+_snn_res.", "", cluster)

    markers=FindAllMarkers(subobj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    markers=markers[markers$p_val_adj < 0.05,]

    write.csv(markers, paste0(cluster_prefix, ".markers.csv"))

    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
    top10genes<-unique(top10$gene)
    
    data.norm=get_seurat_average_expression(subobj, cluster)
    predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
    layer_ids<-names(predict_celltype$max_cta)
    names(layer_ids) <- colnames(data.norm)
    
    seurat_cur_layer<-"seurat_celltype"
    cur_layer="cur_layer"

    oldcluster<-unlist(subobj[[cluster]][[1]])
    newct<-layer_ids[oldcluster]
    cts<-data.frame("seurat_clusters"=oldcluster, cur_layer=newct)
    cts[,seurat_cur_layer] = paste0(cts$seurat_clusters, ": ", cts[,cur_layer])

    ct<-cts[,c("seurat_clusters", seurat_cur_layer)]
    ct<-unique(ct)
    ct<-ct[order(ct$seurat_clusters),]
    
    subobj<-AddMetaData(subobj, cts$seurat_clusters, col.name = "seurat_clusters")
    subobj<-AddMetaData(subobj, cts[,cur_layer], col.name = cur_layer)
    subobj<-AddMetaData(subobj, factor(cts[,seurat_cur_layer], levels=ct[,seurat_cur_layer]), col.name = seurat_cur_layer)
    
    g<-DimPlot(subobj, group.by = "seurat_clusters", label=T) + ggtitle(paste0(pct, ": res", cur_resolution) ) + scale_color_discrete(labels = ct[,seurat_cur_layer])
    if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
      layout <- "ABB"
      g2<-get_bubble_plot(subobj, "seurat_clusters", cur_layer, bubblemap_file, assay="RNA")
      g<-g+g2+plot_layout(design=layout)
      width=6300
    }else{
      width=2300
    }
    
    png(paste0(cluster_prefix, ".umap.png"), width=width, height=2000, res=300)
    print(g)
    dev.off()
    
    width<-max(3000, min(10000, length(unique(subobj$seurat_clusters)) * 150 + 1000))
    height<-max(3000, min(10000, length(top10genes) * 60 + 1000))
    
    subobj<-myScaleData(subobj, top10genes, "RNA")
    g<-DoHeatmap(subobj, assay="RNA", features = top10genes, group.by = seurat_cur_layer, angle = 90) + NoLegend()
    png(paste0(cluster_prefix, ".top10.heatmap.png"), width=3000, height=3000, res=300)
    print(g)
    dev.off()
  }
}
