rm(list=ls()) 
outFile='P10940'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parFile1='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge/result/P10940.final.rds'
parFile2='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge_dr0.2_1_call/result/P10940.scDynamic.meta.rds'
parFile3='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/essential_genes/result/P10940.txt'


setwd('/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge_dr0.2_2_subcluster_rh/result')

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
has_bubblemap <- !is_file_empty(bubblemap_file)
if(has_bubblemap){
  cat("Using bubblemap file ", bubblemap_file, "\n")
}else{
  cat("No bubblemap file\n")
}

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

meta<-readRDS(parFile2)

bHasSignacX<-FALSE
if(exists("parSampleFile4")){
  meta = fill_meta_info_list(parSampleFile4, meta, "signacx_CellStates", "SignacX")
  bHasSignacX<-TRUE
}

bHasSingleR<-FALSE
if(exists('parSampleFile5')){
  meta = fill_meta_info_list(parSampleFile5, meta, "SingleR_labels", "SingleR")
  bHasSingleR<-TRUE
}

obj<-read_object(parFile1)
stopifnot(all(colnames(obj) == rownames(meta)))
obj@meta.data = meta
Idents(obj)<-previous_layer

antibody_bubblemap_file=myoptions$antibody_bubblemap_file
has_antibody_bubblemap <- FALSE
if("ADT" %in% names(obj)){
  has_antibody_bubblemap <- !is_file_empty(antibody_bubblemap_file)
  if(has_antibody_bubblemap){
    cat("Using antibody bubblemap file ", antibody_bubblemap_file, "\n")
  }else{
    antibody_bubblemap_file = NULL
    cat("No antibody bubblemap file\n")
  }
  writeLines(colnames(obj[["ADT"]]), paste0(outFile, ".antibody_markers.txt"))
}

if(by_harmony){
  if(!('batch' %in% colnames(obj@meta.data))){
    obj$batch<-obj$orig.ident
  }
}

if(!is_file_empty(parSampleFile2)){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  obj<-obj[!(rownames(obj) %in% ignore_genes),]
}

draw_dim_plot<-function(obj, previous_layer, file_path){
  old_cts = unique(unlist(obj[[previous_layer]]))
  if(length(old_cts) > 12){
    ncol = 2
  }else{
    ncol = 1
  }
  max_char = max(unlist(lapply(as.character(old_cts), nchar))) + 2
  width = 2000 + max_char * 25 * ncol

  g<-get_dim_plot_labelby(obj, label.by=previous_layer, reduction="umap", label.size = 4, legend.title="") + 
        theme(legend.text = element_text(size = 10)) + guides(fill=guide_legend(ncol=ncol))
  ggsave(file_path, g, width=width, height=2000, dpi=300, units="px", bg="white")
  rm(g)
}

if(!is_file_empty(parSampleFile3)){
  draw_dim_plot(obj, previous_layer, paste0(outFile, ".pre_rename.umap.png"))

  rename_map = read.table(parSampleFile3, sep="\t", header=F)

  meta[,previous_layer]<-as.character(meta[,previous_layer])

  keys = unique(rename_map$V3)
  if("from" %in% rename_map$V2){
    rname = keys[1]
    for(rname in keys){
      rmap = rename_map[rename_map$V3 == rname,]
      from = rmap$V1[rmap$V2=="from"]
      if(!(from %in% unlist(meta[,previous_layer]))){
        stop(paste0("Cannot find ", from, " in obj cell type layer ", previous_layer))
      }
      submeta<-meta[meta[,previous_layer] == from,]

      cluster = rmap$V1[rmap$V2=="cluster"]
      to = rmap$V1[rmap$V2=="to"]
      cluster_column = ifelse('column' %in% rmap$V2, rmap$V1[rmap$V2=="column"], previous_cluster)

      if(all(cluster == "-1")){
        cells<-rownames(submeta)
      }else{
        cur_custers = unique(submeta[,cluster_column])
        if(!(all(cluster %in% unlist(submeta[,cluster_column])))){
          stop(paste0("Cannot find cluster ", paste0(cluster, collapse = "/"), " in cell type ", from, " of cluster ", cluster_column))
        }

        cells<-rownames(submeta)[submeta[,cluster_column] %in% cluster]
      }
      meta[cells,previous_layer]<-to
    }
  }else{
    rname = keys[1]
    for(rname in keys){
      rmap = rename_map[rename_map$V3 == rname,]
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
  
  if("DELETE" %in% meta[,previous_layer]){
    cells = rownames(meta)[meta[,previous_layer] != "DELETE"]
    obj = subset(obj, cells=cells)
  }
  meta=obj@meta.data
  meta[,previous_layer]<-factor_by_count(meta[,previous_layer])
  obj@meta.data = meta

  draw_dim_plot(obj, previous_layer, paste0(outFile, ".post_rename.umap.png"))
}

saveRDS(meta, paste0(outFile, ".meta.rds"))

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes, species = myoptions$species)
  essential_genes<-unique(c(essential_genes, genes_df$gene))

  if(file.exists('fileList6.txt')){
    bubble_file_map = read_file_map('fileList6.txt')
    for(bffile in bubble_file_map){
      if(!file.exists(bffile)){
        stop(paste0("Cannot find bubble file ", bffile))
      }
      genes_df <- read_bubble_genes(bffile, allgenes, species = myoptions$species)
      essential_genes<-unique(c(essential_genes, genes_df$gene))
    }
  }else{
    bubble_file_map = c()
  }
}

default_resolutions=c(seq(from = 0.01, to = 0.09, by = 0.01), seq(from = 0.1, to = 0.5, by = 0.05))
if(is.null(myoptions$resolutions)){
  resolutions=default_resolutions
}else if (myoptions$resolutions == "") {
  resolutions=default_resolutions
}else{
  resolutions=as.numeric(unlist(strsplit(myoptions$resolutions, ",")))
}

layer4map<-split(tiers$Layer4, tiers$Celltype.name)

meta = obj@meta.data
curprefix = prefix

clusters_prefix = paste0(assay, "_snn_res.", format(resolutions, nsmall = 2))
clusters=paste0(assay, "_snn_res.", resolutions)
names(clusters_prefix) = clusters

meta[[previous_layer]]<-factor_by_count(meta[[previous_layer]])
previous_celltypes<-levels(meta[[previous_layer]])
writeLines(previous_celltypes, paste0(outFile, ".cell_types.txt"))
#previous_celltypes<-c("B cells")

DefaultAssay(obj)<-assay

obj$dot<-obj[[previous_layer]]
ctnames<-unique(obj$dot)
ctmap=c(1:length(ctnames))
names(ctmap)<-ctnames
obj@meta.data$dot_cluster<-unlist(ctmap[obj$dot]) + 1000

cur_folder = getwd()
tmp_folder = paste0(cur_folder, "/details")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}
setwd(tmp_folder)

get_bubblemap_file<-function(pct, bubble_file_map, bubblemap_file){
  if(pct %in% names(bubble_file_map)){
    return(bubble_file_map[[pct]])
  }else{
    return(bubblemap_file)
  }
}

filelist<-NULL
allmarkers<-NULL
allcts<-NULL
cluster_index=0
pct<-previous_celltypes[1]
cat("memory used: ", lobstr_mem_used(), "\n")
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
  
  curprefix<-paste0(prefix, ".", celltype_to_filename(pct))

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
                        reduction.name = subumap)
  
  reductions_rds = paste0(curprefix, ".reductions.rds")
  saveRDS(subobj@reductions, reductions_rds)

  min_percentage = 0.05
  bHasCurrentSignacX<-FALSE
  if(bHasSignacX){
    sx<-table(subobj$SignacX)
    sx<-sx[sx > max(5, ncol(subobj) * min_percentage)]
    sxnames<-names(sx)
    
    sxobj<-subset(subobj, SignacX %in% sxnames)
    sxobj$SignacX<-as.character(sxobj$SignacX)

    g4<-get_dim_plot_labelby(sxobj, reduction="umap", label.by = "SignacX", label=T) + ggtitle("SignacX")
    g5<-get_dim_plot_labelby(sxobj, reduction=subumap, label.by = "SignacX", label=T) + ggtitle(paste0(subumap, ": SignacX"))
    g<-g4+g5+plot_layout(ncol=2)

    signacx_file=paste0(curprefix, ".SignacX.umap.png")
    ggsave(signacx_file, g, width=3600, height=1500, dpi=300, units="px", bg="white")
    rm(g4, g5, g)

    bHasCurrentSignacX = any(sxnames != "Unclassified")

    cur_df = data.frame("file"=signacx_file, "type"="SignacX", "resolution"="", "celltype"=pct)
    filelist<-rbind(filelist, cur_df)

    if(has_bubblemap){
      g<-get_sub_bubble_plot(obj, "SignacX", sxobj, "SignacX", bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(curprefix, ".dot.SignacX.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(sxobj, "SignacX"), dpi=300, units="px", bg="white")
      rm(g)

      filelist<-rbind(filelist, c(paste0(getwd(), "/", dot_file), "SignacX celltype", "", pct))
    }

    if(pct %in% names(bubble_file_map)){
      cur_bubblemap_file = bubble_file_map[[pct]]
      g<-get_sub_bubble_plot(obj, "SignacX", sxobj, "SignacX", cur_bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(curprefix, ".dot.SignacX.subcelltypes.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(sxobj, "SignacX"), dpi=300, units="px", bg="white")
      rm(g)

      filelist<-rbind(filelist, c(paste0(getwd(), "/", dot_file), "SignacX sub celltype", "", pct))
    }

    rm(sxobj)
  }

  bHasCurrentSingleR<-FALSE
  if(bHasSingleR){
    sx<-table(subobj$SingleR)
    sx<-sx[sx > max(5, ncol(subobj) * min_percentage)]
    sxnames<-names(sx)

    srobj<-subset(subobj, SingleR %in% sxnames)
    srobj$SingleR<-as.character(srobj$SingleR)
    bHasCurrentSingleR<-TRUE

    g4<-get_dim_plot_labelby(srobj, reduction="umap", label.by = "SingleR", label=T) + ggtitle("SingleR")
    g5<-get_dim_plot_labelby(srobj, reduction=subumap, label.by = "SingleR", label=T) + ggtitle(paste0(subumap, ": SingleR"))
    g<-g4+g5+plot_layout(ncol=2)
    signacr_file=paste0(curprefix, ".SingleR.umap.png")
    ggsave(signacr_file, g, width=3600, height=1500, dpi=300, units="px", bg="white")

    rm(g4, g5, g)

    cur_df = data.frame("file"=signacr_file, "type"="SingleR", "resolution"="", "celltype"=pct)
    filelist<-rbind(filelist, cur_df)

    if(has_bubblemap){
      g<-get_sub_bubble_plot(obj, "SingleR", srobj, "SingleR", bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(curprefix, ".dot.SingleR.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(srobj, "SingleR"), dpi=300, units="px", bg="white")
      rm(g)

      filelist<-rbind(filelist, c(paste0(getwd(), "/", dot_file), "SingleR celltype", "", pct))
    }

    if(pct %in% names(bubble_file_map)){
      cur_bubblemap_file = bubble_file_map[[pct]]
      g<-get_sub_bubble_plot(obj, "SingleR", srobj, "SingleR", cur_bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(curprefix, ".dot.SingleR.subcelltypes.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(srobj, "SingleR"), dpi=300, units="px", bg="white")
      rm(g)

      filelist<-rbind(filelist, c(paste0(getwd(), "/", dot_file), "SingleR sub celltype", "", pct))
    }

    bHasCurrentSingleR = any(sxnames != "unclassified")

    rm(srobj)
  }
  
  cluster = clusters[10]
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
    
    cat(key, "Find marker genes\n")
    markers=FindAllMarkers(subobj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
    markers=markers[markers$p_val_adj < 0.05,]
    
    if(nrow(markers) == 0){
      cat("    no DE gene, pass\n")
      next
    }
    
    markers_file=paste0(cluster_prefix, ".markers.csv")
    write.csv(markers, markers_file)
    
    top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
    top10genes<-unique(top10$gene)
    
    #using global normalized data for cell type annotation
    subobj2<-subset(obj, cells=cells)
    subobj2@meta.data <- subobj@meta.data

    cat(key, "Cell type annotation\n")
    data.norm=get_seurat_average_expression(subobj2, cluster)

    predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
    saveRDS(predict_celltype, paste0(cluster_prefix, ".cta.rds"))
    if(length(predict_celltype$max_cta) > 1){
      Plot_predictcelltype_ggplot2( predict_celltype, 
                            filename=paste0(cluster_prefix, ".cta.png"))
    }

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
    
    # if(bHasCurrentSignacX){
    #   sx<-table(subobj$SignacX)
    #   sx<-sx[sx > max(5, ncol(subobj) * 0.01)]
    #   sxnames<-names(sx)
    #   sxnames<-sxnames[sxnames != "Unclassified"]
    #   sxobj<-subset(subobj, SignacX %in% sxnames)
    #   sxobj$SignacX<-as.character(sxobj$SignacX)
    # }

    # if(bHasCurrentSingleR){
    #   sx<-table(subobj$SingleR)
    #   sx<-sx[sx > max(5, ncol(subobj) * 0.01)]
    #   sxnames<-names(sx)
    #   sxnames<-sxnames[sxnames != "unclassified"]
    #   srobj<-subset(subobj, SingleR %in% sxnames)
    #   srobj$SingleR<-as.character(srobj$SingleR)
    # }

    bar_file=paste0(cluster_prefix, ".bar.png")
    g<-get_barplot(
      ct_meta=subobj@meta.data, 
      bar_file=bar_file,
      cluster_name="display_layer", 
      validation_columns=c("orig.ident", "SignacX", "SingleR"),
      calc_height_per_cluster=250, 
      calc_width_per_cell=50)
    rm(g)

    cat("draw umap file\n")
    g0<-MyDimPlot(obj, label=F, cells.highlight=cells, order = TRUE) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))
    g1<-MyDimPlot(subobj, reduction=subumap, group.by = "orig.ident", label=F) + ggtitle(paste0(pct, ": sample"))
    if(length(unique(subobj$orig.ident)) > 10){
      g1<-g1+NoLegend()
    }

    g2<-MyDimPlot(subobj, reduction="umap", group.by = "seurat_clusters", label=T) + ggtitle(paste0("old umap: res", cur_resolution)) 
    g3<-MyDimPlot(subobj, reduction=subumap, group.by = "seurat_clusters", label=T) + ggtitle(paste0("new umap: res", cur_resolution))
    if(length(unique(subobj$seurat_clusters)) < 20){
      g2<-g2 + scale_color_discrete(labels = ct$display_layer)
      g3<-g3 + scale_color_discrete(labels = ct$display_layer)
    }
    g<-g0+g1+g2+g3+plot_layout(ncol=2)
    umap_file = paste0(cluster_prefix, ".umap.png")
    ggsave(umap_file, g, width=3600, height=3000, dpi=300, units="px", bg="white")
    rm(g, g0, g1, g2, g3)

    cat("draw marker gene heatmap file\n")
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
    ggsave(heatmap_file, gh, width=width, height=height, dpi=300, units="px", bg="white")
    rm(gh)
    
    cur_df = data.frame("file"=paste0(getwd(), "/", c(markers_file, meta_rds, bar_file, umap_file, heatmap_file, reductions_rds)), "type"=c("markers", "meta", "bar", "umap", "heatmap", "reductions"), "resolution"=cur_resolution, "celltype"=pct)

    if(has_bubblemap){
      g<-get_sub_bubble_plot(obj, "dot", subobj, "seurat_celltype", bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(cluster_prefix, ".dot.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(subobj, "seurat_celltype"), dpi=300, units="px", bg="white")
      rm(g)

      cur_df<-rbind(cur_df, c(paste0(getwd(), "/", dot_file), "dot", cur_resolution, pct))
    }

    if(pct %in% names(bubble_file_map)){
      cur_bubblemap_file = bubble_file_map[[pct]]
      g<-get_sub_bubble_plot(obj, "dot", subobj, "seurat_celltype", cur_bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(cluster_prefix, ".dot_celltype_specific.png")
      ggsave(dot_file, width=get_dot_width(g), height=get_dot_height(subobj, "seurat_celltype"), dpi=300, units="px", bg="white")
      rm(g)

      cur_df<-rbind(cur_df, c(paste0(getwd(), "/", dot_file), "dot_celltype_specific", cur_resolution, pct))
    }

    if(has_antibody_bubblemap){
      g<-get_sub_bubble_plot(obj, 
        "dot", 
        subobj, 
        "seurat_celltype", 
        antibody_bubblemap_file, 
        assay="ADT",
        add_num_cell=TRUE, 
        species=NULL)

      dot_file = paste0(cluster_prefix, ".antibody.dot.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(subobj, "seurat_celltype"), dpi=300, units="px", bg="white")
      rm(g)

      cur_df<-rbind(cur_df, c(paste0(getwd(), "/", dot_file), "antibody_dot", cur_resolution, pct))
    }

    filelist<-rbind(filelist, cur_df)
  }
  rm(subobj2)
  rm(subobj)
  cat("after", pct, ", memory used:", lobstr_mem_used(), "\n")
}

rm(obj)
cat("final memory used:", lobstr_mem_used(), "\n")

setwd(cur_folder)

write.csv(filelist, paste0(outFile, ".files.csv"))
