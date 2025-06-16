rm(list=ls()) 
outFile='VK13010_mouse_kidney'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parSampleFile5='fileList5.txt'
parSampleFile7='fileList7.txt'
parFile1='/nobackup/h_cqs/shengq2/temp/20250612_VK13010_scRNA_mouse_kidney/cellbender_nd_seurat_fastmnn/result/VK13010_mouse_kidney.final.rds'
parFile2='/nobackup/h_cqs/shengq2/temp/20250612_VK13010_scRNA_mouse_kidney/cellbender_nd_seurat_fastmnn_dr0.1_1_call/result/VK13010_mouse_kidney.scDynamic.meta.rds'
parFile3='/nobackup/h_cqs/shengq2/temp/20250612_VK13010_scRNA_mouse_kidney/essential_genes/result/VK13010_mouse_kidney.txt'


setwd('/nobackup/h_cqs/shengq2/temp/20250612_VK13010_scRNA_mouse_kidney/cellbender_nd_seurat_fastmnn_dr0.1_2_subcluster/result')

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


options(future.globals.maxSize=1024^3*100) #100G
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
redo_fastmnn<-is_one(myoptions$redo_fastmnn, 0)

assay=ifelse(by_sctransform, "SCT", "RNA")

npcs<-as.numeric(myoptions$pca_dims)
thread=as.numeric(myoptions$thread)

species=myoptions$species
markerfile<-myoptions$db_markers_file
annotate_tcell<-is_one(myoptions$annotate_tcell)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)

output_individual_object<-is_one(myoptions$output_individual_object)
save_intermediate_object<-is_one(myoptions$save_intermediate_object)

min_markers<-20

previous_layer<-myoptions$celltype_layer
cat("previous_layer =", previous_layer, "\n")
previous_cluster<-myoptions$celltype_cluster
cat("previous_cluster =", previous_cluster, "\n")

cur_layer<-myoptions$output_layer
cat("output_layer =", cur_layer, "\n")
seurat_cur_layer=paste0("seurat_", cur_layer)
cat("output_layer_cluster =", seurat_cur_layer, "\n")

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

if(file.exists(parFile3)){
  essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1
}else{
  essential_genes=c()
}

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

cat("Reading object", parFile1, "\n")
obj<-read_object(parFile1)

if(parFile2 != ""){
  meta<-readRDS(parFile2)
}else{
  meta<-obj@meta.data
}

bHasSignacX<-FALSE
if("SignacX_column" %in% names(myoptions)){
  meta$SignacX = meta[[myoptions$SignacX_column]]
  bHasSignacX<-TRUE
}else if(exists("parSampleFile4")){
  meta = fill_meta_info_list(parSampleFile4, meta, "signacx_CellStates", "SignacX")
  bHasSignacX<-TRUE
}else if("SignacX" %in% colnames(obj@meta.data)){
  bHasSignacX<-TRUE
}

bHasSingleR<-FALSE
if("SingleR_column" %in% names(myoptions)){
  meta$SingleR = meta[[myoptions$SingleR_column]]
  bHasSingleR<-TRUE
}else if(exists('parSampleFile5')){
  meta = fill_meta_info_list(parSampleFile5, meta, "SingleR_labels", "SingleR")
  bHasSingleR<-TRUE
}else if("SingleR" %in% colnames(obj@meta.data)){
  bHasSingleR<-TRUE
}

bHasAzimuth<-FALSE
if("Azimuth_column" %in% names(myoptions)){
  meta$Azimuth = meta[[myoptions$Azimuth_column]]
  bHasAzimuth<-TRUE
}else if(exists('parSampleFile7')){
  meta = fill_meta_info_list(parSampleFile7, meta, "Azimuth_finest", "Azimuth")
  bHasAzimuth<-TRUE
}else if("Azimuth" %in% colnames(obj@meta.data)){
  bHasAzimuth<-TRUE
}

stopifnot(all(colnames(obj) == rownames(meta)))
obj@meta.data = meta

if(!previous_layer %in% colnames(obj@meta.data)){
  stop(paste0("Cannot find cell type layer ", previous_layer, " in object meta data. Please check your parameter setting."))
}

if(!previous_cluster %in% colnames(obj@meta.data)){
  stop(paste0("Cannot find cell type layer cluster ", previous_cluster, " in object meta data. Please check your parameter setting."))
}
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

if(!('batch' %in% colnames(obj@meta.data))){
  obj$batch<-obj$orig.ident
}

if(!is_file_empty(parSampleFile2)){
  ignore_gene_files=read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  ignore_genes=unlist(lapply(ignore_gene_files$V1, function(x){
    readLines(x)
  }))
  cat("removing", length(ignore_genes), "genes in file", parSampleFile2, " ...\n")
  #We want to ignore those genes in subclustering but not remove them from whole object.
  #Since we will not save the object in subcluster task, removing those genes here would be most easy way.
  #In subcluster choose task, it would load original object with those genes.
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

meta = obj@meta.data
if(!is_file_empty(parSampleFile3)){
  draw_dim_plot(obj, previous_layer, paste0(outFile, ".pre_rename.umap.png"))

  rename_map = read.table(parSampleFile3, sep="\t", header=F)

  old_layer = previous_layer
  previous_layer=paste0(old_layer, "_renamed")

  meta[,previous_layer]<-as.character(meta[,old_layer])
  if(any(is.na(meta[,previous_layer]))){
    meta[is.na(meta[,previous_layer]),previous_layer] = "DELETE"
  }

  meta$cur_layer = meta[[previous_layer]]
  meta$seurat_clusters = meta[[previous_cluster]]

  all_ct_tbl=rename_map |>
    dplyr::filter(V2 == "ACTIONS")
 
  if(nrow(all_ct_tbl) > 0) {
    all_cts = unique(all_ct_tbl$V3)
    missed_cts = setdiff(all_cts, unique(unlist(meta[,previous_layer])))
    if(length(missed_cts) > 0){
      current_cts = unique(unlist(meta[,previous_layer]))
      stop(paste0("Cannot find cell types ", paste0(missed_cts, collapse = "/"), " in obj cell type layer ", previous_layer, ": ", paste0(current_cts, collapse = "/")))
    }

    move_ct_tbl=all_ct_tbl |> dplyr::filter(grepl("^MOVE", V1))
    if(nrow(move_ct_tbl) > 0) {
      ct=move_ct_tbl$V3[1]
      move_cts=sort(unique(move_ct_tbl$V3))
      for(ct in move_cts){
        print(ct)
        ct_tbl=move_ct_tbl |> dplyr::filter(V3 == ct)
        cur_meta = meta[meta$cur_layer == ct,]
        cur_meta = process_actions(ct_tbl, cur_meta)
        meta[rownames(cur_meta),previous_layer] = cur_meta$cur_layer
      }
    }

    other_ct_tbl=all_ct_tbl |> dplyr::filter(!grepl("^MOVE", V1))
    if(nrow(other_ct_tbl) > 0){
      ct=other_ct_tbl$V3[1]
      other_cts=sort(unique(other_ct_tbl$V3))
      for(ct in other_cts){
        print(ct)
        ct_tbl=other_ct_tbl |> dplyr::filter(V3 == ct)
        cur_meta = meta[meta$cur_layer == ct,]
        cur_meta = process_actions(ct_tbl, cur_meta)
        cur_meta[cur_meta$seurat_clusters < 0, "cur_layer"] = "DELETE"
        meta[rownames(cur_meta),previous_layer] = cur_meta$cur_layer
      }
    }
    
    meta = meta |> dplyr::select(-cur_layer, -seurat_clusters)
  }

  rename_map=rename_map |>
    dplyr::filter(V2 != "ACTIONS")

  keys = unique(rename_map$V3)
  if("from" %in% rename_map$V2){
    rname = keys[1]
    for(rname in keys){
      cat("renaming", rname, "\n")
      rmap = rename_map[rename_map$V3 == rname,]
      from = rmap$V1[rmap$V2=="from"]
      if(!(from %in% unique(unlist(meta[,previous_layer])))){
        stop(paste0("Cannot find ", from, " in obj cell type layer ", previous_layer))
      }
      submeta<-meta[meta[,previous_layer] == from,]

      cluster = rmap$V1[rmap$V2=="cluster"]
      to = rmap$V1[rmap$V2=="to"]
      
      if(all(cluster == "-1")){
        cells<-rownames(submeta)
      }else{
        cluster_column = ifelse('column' %in% rmap$V2, rmap$V1[rmap$V2=="column"], previous_cluster)
        if(!(cluster_column %in% colnames(submeta))){
          stop(paste0("Cannot find column ", cluster_column, " in rename configuration key=", rname))
        }

        cur_custers = unique(submeta[,cluster_column])
        missed_clusters = setdiff(cluster, cur_custers)
        if(length(missed_clusters) > 0){
          stop(paste0("Cannot find cluster ", paste0(missed_clusters, collapse = "/"), " in cell type ", from, " of cluster ", cluster_column))
        }

        cells<-rownames(submeta)[submeta[,cluster_column] %in% cluster]
      }
      meta[cells,previous_layer]<-to
      cat("renaming", rname, "done\n")
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

  cat("renaming done\n")

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
tmp_folder = paste0(cur_folder, "/details/")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}

get_bubblemap_file<-function(pct, bubble_file_map, bubblemap_file){
  if(pct %in% names(bubble_file_map)){
    return(bubble_file_map[[pct]])
  }else{
    return(bubblemap_file)
  }
}

check_cell_type<-function(subobj, ct_column, filelist, pct, curprefix, species, subumap, has_bubblemap, bubblemap_file, bubble_file_map){
  stopifnot(ct_column %in% colnames(subobj@meta.data))

  sxobj=get_filtered_obj(subobj, ct_column)
  sxnames<-levels(sxobj[[ct_column]])

  g4<-get_dim_plot_labelby(sxobj, reduction="umap", label.by = ct_column, label=T) + ggtitle(ct_column)
  g5<-get_dim_plot_labelby(sxobj, reduction=subumap, label.by = ct_column, label=T) + ggtitle(paste0(subumap, ": ", ct_column))
  g<-g4+g5+plot_layout(ncol=2)

  nlens=nchar(sxobj@meta.data[[ct_column]])
  nlens=nlens[!is.na(nlens)]
  max_len=max(nlens)
  width=ifelse(max_len < 30, 3600, 4800)

  ct_file=paste0(curprefix, ".", ct_column, ".umap.png")
  ggsave(ct_file, g, width=width, height=1500, dpi=300, units="px", bg="white")
  rm(g4, g5, g)

  cur_df = data.frame("file"=ct_file, "type"=ct_column, "resolution"="", "celltype"=pct)
  filelist<-rbind(filelist, cur_df)

  if(has_bubblemap){
    g<-get_sub_bubble_plot(obj, ct_column, sxobj, ct_column, bubblemap_file, add_num_cell=TRUE, species=species)

    dot_file = paste0(curprefix, ".dot.", ct_column, ".png")
    ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(sxobj, ct_column), dpi=300, units="px", bg="white")
    rm(g)

    filelist<-rbind(filelist, c(dot_file, paste0(ct_column, " celltype"), "", pct))
  }

  if(pct %in% names(bubble_file_map)){
    cur_bubblemap_file = bubble_file_map[[pct]]
    g<-get_sub_bubble_plot(obj, ct_column, sxobj, ct_column, cur_bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

    dot_file = paste0(curprefix, ".dot.", ct_column, ".subcelltypes.png")
    ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(sxobj, ct_column), dpi=300, units="px", bg="white")
    rm(g)

    filelist<-rbind(filelist, c(dot_file, paste0(ct_column, " sub celltype"), "", pct))
  }

  rm(sxobj)

  return(filelist)
}

filelist<-NULL
allmarkers<-NULL
allcts<-NULL
cluster_index=0
#previous_celltypes=rev(previous_celltypes)
pct<-previous_celltypes[1]
cat("memory used: ", lobstr_mem_used(), "\n")
for(pct in previous_celltypes){
  key = paste0(previous_layer, ": ", pct, ":")
  cat(key, "Start\n")
  cells<-rownames(meta)[meta[,previous_layer] == pct]
  subobj<-subset(obj, cells=cells)
  
  stopifnot(all(subobj[[previous_layer]] == pct))
  
  pca_npcs<-min(round(length(cells)/2), 50)
  cur_npcs=min(pca_npcs, npcs)
  cur_pca_dims=1:cur_npcs

  curreduction = reduction

  k_n_neighbors<-min(cur_npcs, 20)
  u_n_neighbors<-min(cur_npcs, 30)
  
  curprefix<-paste0(tmp_folder, prefix, ".", celltype_to_filename(pct))

  subumap = "subumap"
  subobj = sub_cluster( subobj, 
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
                        redo_fastmnn = redo_fastmnn,
                        thread=thread,
                        detail_prefix=curprefix)
  
  cat("saving reductions ...\n")
  reductions_rds = paste0(curprefix, ".reductions.rds")
  saveRDS(subobj@reductions, reductions_rds)

  validation_columns=c("orig.ident")
  if(bHasSignacX){
    filelist = check_cell_type(subobj, "SignacX", filelist, pct, curprefix, species, subumap, has_bubblemap, bubblemap_file, bubble_file_map)
    validation_columns = c(validation_columns, "SignacX")
  }

  if(bHasSingleR){
    filelist = check_cell_type(subobj, "SingleR", filelist, pct, curprefix, species, subumap, has_bubblemap, bubblemap_file, bubble_file_map)
    validation_columns = c(validation_columns, "SingleR")
  }

  if(bHasAzimuth){
    filelist = check_cell_type(subobj, "Azimuth", filelist, pct, curprefix, species, subumap, has_bubblemap, bubblemap_file, bubble_file_map)
    validation_columns = c(validation_columns, "Azimuth")
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
    
    if(nrow(markers) < 2){
      #with less than 2 DE genes, we cannot draw heatmap.
      cat("    less than 2 DE gene, pass\n")
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

    bar_file=paste0(cluster_prefix, ".bar.png")
    g<-get_barplot(
      ct_meta=subobj@meta.data, 
      bar_file=bar_file,
      cluster_name="display_layer", 
      validation_columns=validation_columns,
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
    
    cur_df = data.frame("file"=c(markers_file, meta_rds, bar_file, umap_file, heatmap_file, reductions_rds), "type"=c("markers", "meta", "bar", "umap", "heatmap", "reductions"), "resolution"=cur_resolution, "celltype"=pct)

    if(has_bubblemap){
      g<-get_sub_bubble_plot(obj, "dot", subobj, "seurat_celltype", bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(cluster_prefix, ".dot.png")
      ggsave(dot_file, g, width=get_dot_width(g), height=get_dot_height(subobj, "seurat_celltype"), dpi=300, units="px", bg="white")
      rm(g)

      cur_df<-rbind(cur_df, c(dot_file, "dot", cur_resolution, pct))
    }

    if(pct %in% names(bubble_file_map)){
      cur_bubblemap_file = bubble_file_map[[pct]]
      g<-get_sub_bubble_plot(obj, "dot", subobj, "seurat_celltype", cur_bubblemap_file, add_num_cell=TRUE, species=myoptions$species)

      dot_file = paste0(cluster_prefix, ".dot_celltype_specific.png")
      ggsave(dot_file, width=get_dot_width(g), height=get_dot_height(subobj, "seurat_celltype"), dpi=300, units="px", bg="white")
      rm(g)

      cur_df<-rbind(cur_df, c(dot_file, "dot_celltype_specific", cur_resolution, pct))
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

      cur_df<-rbind(cur_df, c(dot_file, "antibody_dot", cur_resolution, pct))
    }

    filelist<-rbind(filelist, cur_df)
  }
  rm(subobj2)
  cat("after", pct, ", memory used:", lobstr_mem_used(), "\n")

  if(output_individual_object){
    saveRDS(subobj, paste0(curprefix, ".", pct, ".obj.rds"))
  }
  rm(subobj)
}

rm(obj)
cat("final memory used:", lobstr_mem_used(), "\n")

write.csv(filelist, paste0(outFile, ".files.csv"))
