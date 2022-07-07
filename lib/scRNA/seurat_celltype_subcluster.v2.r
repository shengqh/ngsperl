rm(list=ls()) 
outFile='mouse_8363'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/scratch/jbrown_lab/shengq2/projects/20220630_scRNA_8363_mouse/seurat_merge/result/mouse_8363.final.rds'
parFile2='C:/projects/scratch/jbrown_lab/shengq2/projects/20220630_scRNA_8363_mouse/seurat_merge_01_dynamic/result/mouse_8363.scDynamic.meta.rds'
parFile3='C:/projects/scratch/jbrown_lab/shengq2/projects/20220630_scRNA_8363_mouse/essential_genes/result/mouse_8363.txt'


setwd('C:/projects/scratch/jbrown_lab/shengq2/projects/20220630_scRNA_8363_mouse/seurat_merge_02_subcluster.v2/result')

### Parameter setting end ###

source("scRNA_func.r")
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
library(testit)
library(stringr)

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

bHasSignacX<-FALSE
if(exists("parFile4")){
  if(parFile4 != ""){
    signacX<-readRDS(parFile4)
    assert(rownames(signacX) == rownames(obj@meta.data))
  
    ct_map=c('T.CD4.memory'='T.CD4', 
      'T.CD4.naive'='T.CD4', 
      'T.CD8.cm'='T.CD8',
      'T.CD8.em'='T.CD8',
      'T.CD8.naive'='T.CD8')
    signacX$signacx_CellStates_slim<-as.character(signacX$signacx_CellStates)
    for(ct_name in names(ct_map)){
      signacX$signacx_CellStates_slim[signacX$signacx_CellStates_slim==ct_name]=ct_map[ct_name]
    }
  
    obj<-AddMetaData(obj, signacX$signacx_CellStates_slim, col.name = "signacx_CellStates")
    bHasSignacX<-TRUE
  }
}

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  bubble_genes<-unique(genes_df$gene)
}

resolutions=c(seq(from = 0.01, to = 0.09, by = 0.01), seq(from = 0.1, to = 0.5, by = 0.1))

layer4map<-split(tiers$Layer4, tiers$Celltype.name)

meta = obj@meta.data
curprefix = prefix

clusters_prefix = paste0(assay, "_snn_res.", format(resolutions, nsmall = 2))
clusters=paste0(assay, "_snn_res.", resolutions)
names(clusters_prefix) = clusters

celltypes<-unlist(meta[[previous_layer]])
tblct<-table(celltypes)
tblct<-tblct[order(tblct, decreasing = F)]

previous_celltypes<-names(tblct)
#previous_celltypes<-c("B cells")

DefaultAssay(obj)<-assay

g00<-DimPlot(obj, group.by=previous_layer, label=T)
g01<-DimPlot(obj, group.by="orig.ident", label=T)

filelist<-NULL
allmarkers<-NULL
allcts<-NULL
cluster_index=0
pct<-previous_celltypes[2]
for(pct in previous_celltypes){
  key = paste0(previous_layer, ": ", pct, ":")
  cells<-rownames(meta)[meta[,previous_layer] == pct]
  subobj<-subset(obj, cells=cells)
  
  stopifnot(all(subobj[[previous_layer]] == pct))
  
  pca_npcs<-min(round(length(cells)/2), 50)
  cur_npcs=min(pca_npcs, npcs)
  cur_pca_dims=1:cur_npcs

  k_n_neighbors<-min(cur_npcs, 20)
  u_n_neighbors<-min(cur_npcs, 30)
  
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
    cat(key, "RunPCA\n")
    subobj<-RunPCA(subobj, npcs=cur_npcs)
    curreduction="pca"
  }

  DefaultAssay(subobj)<-assay

  cat(key, "FindClusters\n")
  subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)
  subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolutions, verbose=FALSE)
  
  
  g11<-DimPlot(subobj, reduction="umap", group.by = "orig.ident", label=F)
  
  cat(key, "RunUMAP\n")
  umap_names<-c()
  ops<-data.frame("nn"=c(30,20,10), "min.dist"=c(0.3,0.1,0.05))
  for(idx in c(1:nrow(ops))){
    nn=ops[idx, 'nn']
    nn<-min(nn, u_n_neighbors)
    min.dist=ops[idx, 'min.dist']
    umap_name = paste0("umap_nn", nn, "_dist", min.dist)
    cat(umap_name, "\n")
    umap_names<-c(umap_names, umap_name)
    umap_key = paste0("UMAPnn", nn, "dist", min.dist * 100, "_")
    subobj<-RunUMAP(object = subobj, reduction=curreduction, reduction.key=umap_key, reduction.name=umap_name, n.neighbors=nn, min.dist=min.dist, dims=cur_pca_dims, verbose = FALSE)
  }
# 
#   g<-NULL
#   for(umap_name in umap_names){  
#     cat(umap_name, "\n")
#     cg<-DimPlot(subobj, reduction = umap_name, group.by="orig.ident") + ggtitle(umap_name)
#     if(is.null(g)){
#       g<-cg
#     }else{
#       g<-g+cg
#     }
#   }
#   g<-g+plot_layout(ncol=3)
#   png(paste0(curprefix, ".iumap.png"), width=6600, height=6000, res=300)
#   print(g)
#   dev.off()
  
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
    
    if(nrow(markers) == 0){
      cat("    only DE gene, pass\n")
      next
    }
    
    markers_file=paste0(cluster_prefix, ".markers.csv")
    write.csv(markers, markers_file)
    
    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
    top10genes<-unique(top10$gene)
    
    #using global normalized data for cell type annotation
    subobj2<-subset(obj, cells=cells)
    subobj2@meta.data <- subobj@meta.data

    data.norm=get_seurat_average_expression(subobj2, cluster)
    predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
    layer_ids<-names(predict_celltype$max_cta)
    names(layer_ids) <- colnames(data.norm)
    
    seurat_cur_layer<-"seurat_celltype"
    cur_layer="cur_layer"
    
    oldcluster<-unlist(subobj[[cluster]][[1]])
    newct<-layer_ids[oldcluster]
    cts<-data.frame("seurat_clusters"=oldcluster, cur_layer=newct)
    cts[,seurat_cur_layer] = paste0(cts$seurat_clusters, ": ", cts[,cur_layer])
    new_layer_ids<-word(cts[,seurat_cur_layer], 1,4, sep=" ")
    new_layer_ids[is.na(new_layer_ids)]=cts[is.na(new_layer_ids),seurat_cur_layer]
    cts$display_layer = new_layer_ids

    ct<-cts[,c("seurat_clusters", seurat_cur_layer, "display_layer")]
    ct<-unique(ct)
    ct<-ct[order(ct$seurat_clusters),]
    
    subobj<-AddMetaData(subobj, cts$seurat_clusters, col.name = "seurat_clusters")
    subobj<-AddMetaData(subobj, cts[,cur_layer], col.name = cur_layer)
    subobj<-AddMetaData(subobj, factor(cts[,seurat_cur_layer], levels=ct[,seurat_cur_layer]), col.name = seurat_cur_layer)
    subobj<-AddMetaData(subobj, factor(cts$display_layer, levels=ct$display_layer), col.name = "display_layer")
    
    meta_rds = paste0(cluster_prefix, ".meta.rds")
    saveRDS(subobj@meta.data, meta_rds)
    
    bHasCurrentSignacX<-FALSE
    if(bHasSignacX){
      sx<-table(subobj$signacx_CellStates)
      sx<-sx[sx>5]
      sxnames<-names(sx)
      sxnames<-sxnames[sxnames != "Unclassified"]
      if(length(sxnames)>0){
        sxobj<-subset(subobj, signacx_CellStates %in% sxnames)
        sxobj$signacx_CellStates<-as.character(sxobj$signacx_CellStates)
        bHasCurrentSignacX<-TRUE
      }
    }

    g0<-DimPlot(obj, label=F, cells.highlight =cells) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))
    g1<-DimPlot(subobj, reduction="umap", group.by = "seurat_clusters", label=T) + ggtitle(paste0(pct, ": old umap: res", cur_resolution)) + scale_color_discrete(labels = ct$display_layer)
    if(bHasCurrentSignacX){
      g2<-DimPlot(sxobj, reduction="umap", group.by = "signacx_CellStates", label=F) + ggtitle(paste0(pct, ": signacx"))
    }else{
      g2<-DimPlot(subobj, reduction="umap", group.by = "orig.ident", label=F) + ggtitle(paste0(pct, ": sample"))
    }
    g<-g0+g1+g2+get_groups_dot(subobj, "seurat_clusters", "orig.ident")
    for(umap_name in umap_names){
      gu<-DimPlot(subobj, reduction=umap_name, group.by = "seurat_clusters", label=T) + ggtitle(paste0("res", cur_resolution, ": ", umap_name)) + scale_color_discrete(labels = ct$display_layer)
      g<-g+gu
    }
    if(bHasCurrentSignacX){
      g<-g+get_groups_dot(sxobj, "seurat_clusters", "signacx_CellStates")
      for(umap_name in umap_names){
        gu<-DimPlot(sxobj, reduction=umap_name, group.by = "signacx_CellStates", label=F) + ggtitle(paste0(umap_name, ": signacX"))
        g<-g+gu
      }
      g<-g+get_groups_dot(sxobj, "signacx_CellStates", "orig.ident")
    }

    width=10000
    if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
      #using global normalized data for bubblemap
      subobj2@meta.data <- subobj@meta.data

      gb<-get_bubble_plot(subobj2, "seurat_clusters", cur_layer, bubblemap_file, assay="RNA") + theme(text = element_text(size=20))
      g<-g+gb
      if(bHasCurrentSignacX){
        layout="ABCD
EFGH
IJKL
MMMM"
        height=6000
      }else{
        layout="ABCD
EFGI
HHHH"
        height=4500
      }
    }else{
      if(bHasCurrentSignacX){
        layout="ABCD
EFGH
IJKL"
        height=4500
      }else{
        layout="ABCD
EFGH"
        height=3000
      }
    }
    gu<-g+plot_layout(design=layout)

    #remove subobj2
    rm(subobj2)
    
    umap_file = paste0(cluster_prefix, ".umap.png")
    png(umap_file, width=width, height=height, res=300)
    print(gu)
    dev.off()
    
    width<-max(3000, min(10000, length(unique(subobj$seurat_clusters)) * 150 + 1000))
    height<-max(3000, min(10000, length(top10genes) * 60 + 1000))
    
    subobj<-myScaleData(subobj, top10genes, "RNA")
    gh<-DoHeatmap(subobj, assay="RNA", features = top10genes, group.by = seurat_cur_layer, angle = 90) + NoLegend()
    
    heatmap_file = paste0(cluster_prefix, ".top10.heatmap.png")
    png(heatmap_file, width=3000, height=3000, res=300)
    print(gh)
    dev.off()

#     layout="NNABCD
# NNEFGH
# NNIJKL
# NNMMMM"
#     gg<-g+gh+plot_layout(design=layout)
#     umap_file = paste0(cluster_prefix, ".umap.png")
#     png(umap_file, width=width+height, height=height, res=300)
#     print(gg)
#     dev.off()
    
    cur_df = data.frame("file"=paste0(getwd(), "/", c(markers_file, meta_rds, umap_file, heatmap_file)), "type"=c("markers", "meta", "umap", "heatmap"), "resolution"=cur_resolution, "celltype"=pct)
    filelist<-rbind(filelist, cur_df)
  }
}

write.csv(filelist, paste0(outFile, ".files.csv"))
