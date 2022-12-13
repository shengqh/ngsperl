rm(list=ls()) 
outFile='P9061'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_merge/result/P9061.final.rds'
parFile2='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_merge_multires_01_call/result/P9061.meta.rds'
parFile3='/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/essential_genes/result/P9061.txt'


setwd('/scratch/vickers_lab/projects/20221201_scRNA_9061_mouse/seurat_sct_merge_multires_02_subcluster/result')

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

by_sctransform<-is_one(myoptions$by_sctransform)
by_integration<-is_one(myoptions$by_integration)
reduction<-myoptions$reduction
by_harmony<-reduction=="harmony"
redo_harmony<-is_one(myoptions$redo_harmony, 0)
assay=get_assay(by_sctransform, by_integration, by_harmony)

npcs<-as.numeric(myoptions$pca_dims)

species=myoptions$species
markerfile<-myoptions$db_markers_file
annotate_tcell<-is_one(myoptions$annotate_tcell)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)

min_markers<-20

previous_layer<-myoptions$celltype_layer
previous_cluster<-myoptions$celltype_cluster

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

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = species,
                             curated_markers_file = myoptions$curated_markers_file,
                             HLA_panglao5_file = HLA_panglao5_file,
                             layer="Layer4",
                             remove_subtype_str = "",
                             combined_celltype_file = NULL)

tiers = ctdef$tiers

cell_activity_database<-ctdef$cell_activity_database

prefix<-outFile

if(!exists("obj")){
  obj<-read_object(parFile1, parFile2)
  Idents(obj)<-previous_layer
}

if(!is_file_empty(parSampleFile2)){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  obj<-obj[!(rownames(obj) %in% ignore_genes),]
}

if(!is_file_empty(parSampleFile3)){
  rename_map = read.table(parSampleFile3, sep="\t", header=F)

  meta<-obj@meta.data

  keys = unique(rename_map$V3)
  if("from" %in% rename_map$V2){
    rname = keys[1]
    for(rname in keys){
      rmap = rename_map[rename_map$V3 == rname]
      from = rmap$V1[rmap$V2=="from"]
      cluster = rmap$V1[rmap$V2=="cluster"]
      to = rmap$V1[rmap$V2=="to"]

      if(!(from %in% unlist(meta[,previous_layer]))){
        stop(paste0("Cannot find ", from, " in obj cell type layer ", previous_layer))
      }
      submeta<-meta[meta[,previous_layer] == from,]

      if(!(cluster %in% unlist(submeta[,previous_cluster]))){
        stop(paste0("Cannot find cluster ", cluster, " in cell type ", from, " of cluster ", previous_cluster))
      }

      cells<-rownames(submeta)[submeta[,previous_cluster]==cluster]
      meta[cells,previous_layer]<-to
    }
  }else{
    rname = keys[1]
    for(rname in keys){
      rmap = rename_map[rename_map$V3 == rname]
      cluster = rmap$V1[rmap$V2=="cluster"]
      to = rmap$V1[rmap$V2=="to"]

      if(!(cluster %in% unlist(meta[,previous_cluster]))){
        stop(paste0("Cannot find cluster ", cluster, " in obj cell type cluster ", previous_cluster))
      }

      cells<-rownames(meta)[meta[,previous_cluster]==cluster]
      meta[cells,previous_layer]<-to
    }
  }

  obj@meta.data<-meta
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
ordered_celltypes<-previous_celltypes[order(previous_celltypes)]
writeLines(ordered_celltypes, paste0(outFile, ".cell_types.txt"))
#previous_celltypes<-c("B cells")

DefaultAssay(obj)<-assay

filelist<-NULL
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

  curreduction = reduction

  k_n_neighbors<-min(cur_npcs, 20)
  u_n_neighbors<-min(cur_npcs, 30)
  
  curprefix<-paste0(prefix, ".", celltype_to_filename=(pct))

  subumap = "subumap"
  subobj = sub_cluster(subobj, 
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
                        key,
                        do_umap = TRUE,
                        reduction.name = subumap,
                        previous_layer = NA)
  
  reductions_rds = paste0(curprefix, ".reductions.rds")
  saveRDS(subobj@reductions, reductions_rds)
  
  cat(key, "Find marker genes\n")
  cluster = clusters[5]
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

    # bar file
    bar_file=paste0(cluster_prefix, ".bar.png")
    gb<-get_groups_dot(subobj, "display_layer", "orig.ident")
    if(bHasCurrentSignacX){
      gb<-gb+get_groups_dot(subobj, "display_layer", "signacx_CellStates") + plot_layout(ncol=1)
    }
    height=ifelse(bHasCurrentSignacX, 2200, 1100)
    png(bar_file, width=3000, height=height, res=300)
    print(gb)
    dev.off()

    # umap file
    g0<-DimPlot(obj, label=F, cells.highlight=cells, order = TRUE) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))
    g1<-DimPlot(subobj, reduction="umap", group.by = "seurat_clusters", label=T) + ggtitle(paste0("old umap: res", cur_resolution)) + scale_color_discrete(labels = ct$display_layer)
    g3<-DimPlot(subobj, reduction=subumap, group.by = "orig.ident", label=F) + ggtitle(paste0(pct, ": sample"))
    g4<-DimPlot(subobj, reduction=subumap, group.by = "seurat_clusters", label=T) + scale_color_discrete(labels = ct$display_layer)
    if(bHasCurrentSignacX){
      g2<-DimPlot(sxobj, reduction="umap", group.by = "signacx_CellStates", label=F) + ggtitle("signacx")
      g5<-DimPlot(sxobj, reduction=subumap, group.by = "signacx_CellStates", label=F) + ggtitle(paste0(subumap, ": signacX"))
      g<-g0+g1+g2+g3+g4+g5
    }else{
      g<-g0+g1+g3+g4
    }
    ncol=ifelse(bHasCurrentSignacX, 3, 2)
    g<-g+plot_layout(ncol=ncol)
    width=ncol * 1800
    height=3000
    umap_file = paste0(cluster_prefix, ".umap.png")
    png(umap_file, width=width, height=height, res=300)
    print(g)
    dev.off()

    # marker gene heatmap
    subobj<-myScaleData(subobj, top10genes, "RNA")
    if(ncol(subobj) > 5000){
      subsampled <- subobj[, sample(colnames(subobj), size=5000, replace=F)]
      gh<-DoHeatmap(subsampled, assay="RNA", features = top10genes, group.by = seurat_cur_layer, angle = 90) + NoLegend()
      rm(subsampled)
    }else{
      gh<-DoHeatmap(subobj, assay="RNA", features = top10genes, group.by = seurat_cur_layer, angle = 90) + NoLegend()
    }
    
    heatmap_file = paste0(cluster_prefix, ".top10.heatmap.png")

    width<-max(3000, min(10000, length(unique(subobj$seurat_clusters)) * 150 + 1000))
    height<-max(3000, min(10000, length(top10genes) * 60 + 1000))
    png(heatmap_file, width=width, height=height, res=300)
    print(gh)
    dev.off()
    
    cur_df = data.frame("file"=paste0(getwd(), "/", c(markers_file, meta_rds, bar_file, umap_file, heatmap_file, reductions_rds)), "type"=c("markers", "meta", "bar", "umap", "heatmap", "reductions"), "resolution"=cur_resolution, "celltype"=pct)

    if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
      #using global normalized data for bubblemap
      subobj2@meta.data <- subobj@meta.data

      ncluster<-length(unique(subobj2$seurat_clusters))
      height=max(2000, ncluster*250 + 1000)

      g<-get_bubble_plot(subobj2, "seurat_clusters", cur_layer, bubblemap_file, assay="RNA") + theme(text = element_text(size=20))
      dot_file = paste0(cluster_prefix, ".dot.png")
      png(dot_file, width=6000, height=height, res=300)
      print(g)
      dev.off()

      cur_df<-rbind(cur_df, c(paste0(getwd(), "/", dot_file), "dot", cur_resolution, pct))
    }

    filelist<-rbind(filelist, cur_df)
  }
}

write.csv(filelist, paste0(outFile, ".files.csv"))

library('rmarkdown')
rmarkdown::render("seurat_celltype_subcluster.v3.rmd",output_file=paste0(outFile,".subcluster.html"))

