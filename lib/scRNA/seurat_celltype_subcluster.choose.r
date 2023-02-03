rm(list=ls()) 
outFile='mouse_8870'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1='/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/seurat_sct_harmony/result/mouse_8870.final.rds'
parFile2='/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/seurat_sct_harmony_dr0.2_nrh_01_call/result/mouse_8870.scDynamic.meta.rds'
parFile3='/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/essential_genes/result/mouse_8870.txt'
parFile4='/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/seurat_sct_harmony_dr0.2_nrh_02_subcluster/result/mouse_8870.files.csv'


setwd('/scratch/jbrown_lab/shengq2/projects/20221117_scRNA_8870_mouse/seurat_sct_harmony_dr0.2_nrh_03_choose/result')

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
library(testit)

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

output_heatmap=FALSE

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
reduction<-myoptions$reduction
assay=ifelse(by_sctransform, "SCT", "RNA")

previous_layer<-myoptions$celltype_layer
cur_layer<-myoptions$output_layer
seurat_clusters = "seurat_clusters"
seurat_cur_layer=paste0("seurat_", cur_layer)
resolution_col = "resolution"

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1
essential_genes<-""

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)
bubble_width <-ifelse(is.null(myoptions$bubble_width), 6000, as.numeric(myoptions$bubble_width))

prefix<-outFile

if(!exists("obj")){
  obj<-read_object(parFile1, parFile2)
  obj<-factorize_layer(obj, previous_layer)
  Idents(obj)<-previous_layer
}
# obj<-readRDS("crs.final.rds")

pre_choose_file=paste0(prefix, ".pre.umap.png")

post_rename_umap = file.path(dirname(parFile4), outFile, ".post_rename.umap.png")
if(file.exists(post_rename_umap)){
  file.copy(post_rename_umap, pre_choose_file)
}else{
  obj<-build_dummy_cluster(obj, label.by=previous_layer, new_cluster_name = "old_clusters")
  g<-get_dim_plot(obj, group.by ="old_clusters", label.by="old_clusters_label", label.size = 8, legend.title="") + 
    theme(legend.text = element_text(size = 20)) + ggtitle("")
  png(pre_choose_file, width=3200, height=2000, res=300)
  print(g)
  dev.off()
  obj$old_clusters<-NULL
  obj$old_clusters_label<-NULL
}

obj<-AddMetaData(obj, obj[[previous_layer]], col.name = cur_layer)
obj<-unfactorize_layer(obj, cur_layer)
obj<-AddMetaData(obj, -1, col.name = seurat_clusters)
obj<-AddMetaData(obj, -1, col.name = resolution_col)
obj<-AddMetaData(obj, 0, col.name = "sub_clusters")

#set the subuamp using default umap
#it would be replaced cell type by cell type using real sub umap in subclustering
subumap<-Embeddings(obj, reduction = "umap")
colnames(subumap)<-c("SUBUMAP_1", "SUBUMAP_2")
obj[["subumap"]] <- CreateDimReducObject(embeddings = subumap, key = "SUBUMAP_", assay = DefaultAssay(obj))

clcounts<-table(obj[[cur_layer]])

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
  bubble_genes<-unique(genes_df$gene)
}

best_res_tbl<-read.table(parSampleFile3, sep="\t", header=F, stringsAsFactors = F)

res_files<-read.csv(parFile4, header=T)

#res_files$file<-gsub("/gpfs52", "c:/projects", res_files$file)

if(!("type" %in% colnames(res_files))){
  res_files$type<-unlist(apply(res_files, 1, function(x){
    if(grepl('heatmap', x['file'])){
      return('heatmap')
    }
    if(grepl('markers', x['file'])){
      return('markers')
    }
    if(grepl('umap', x['file'])){
      return('umap')
    }
    if(grepl('meta', x['file'])){
      return('meta')
    }
    stop(x)
  }))
}

if(!all(best_res_tbl$V3 %in% names(clcounts))){
  defined_ct<-unique(best_res_tbl$V3[!(best_res_tbl$V3 %in% names(clcounts))])
  stop(paste0("those cell types were defined at celltype_subclusters_table but not exists:", paste0(defined_ct, sep=",")))
}

if(!all(names(clcounts) %in% best_res_tbl$V3)){
  miss_ct<-unique(names(clcounts)[!(names(clcounts) %in% best_res_tbl$V3)])
  stop(paste0("those cell types were not defined at celltype_subclusters_table :", paste0(miss_ct, collapse=",")))
}

#remove cell type first
remove_cts<-best_res_tbl$V3[best_res_tbl$V2=="resolution" & best_res_tbl$V1 == "-1"]
if(length(remove_cts) > 0){
  cells<-colnames(obj)[!(unlist(obj[[cur_layer]]) %in% remove_cts)]
  obj<-subset(obj, cells=cells)
}

if(output_heatmap){
  #find markers for cell types
  cat("FindAllMarkers ...\n")
  ct_markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
  write.csv(ct_markers, paste0(outFile, ".markers.csv"))
  ct_top10<-get_top10_markers(ct_markers)
  ct_top10_map<-split(ct_top10$gene, ct_top10$cluster)
  cat("FindAllMarkers done ...\n")
}

meta = obj@meta.data
meta[, cur_layer] = as.character(meta[, cur_layer])

curprefix = prefix

celltypes<-unlist(meta[[previous_layer]])
tblct<-table(as.character(celltypes))
tblct<-tblct[order(tblct, decreasing = T)]

previous_celltypes<-names(tblct)
writeLines(previous_celltypes, paste0(outFile, ".orig_cell_types.txt"))
#previous_celltypes<-c("B cells")

DefaultAssay(obj)<-assay

if(output_heatmap){
  allmarkers<-NULL
}

sub_dir<-dirname(parFile4)

cur_folder = getwd()
tmp_folder = paste0(cur_folder, "/details")
if(!dir.exists(tmp_folder)){
  dir.create(tmp_folder)
}
setwd(tmp_folder)

cluster_index=0
pct<-previous_celltypes[1]
for(pct in previous_celltypes){
  cat(pct, "\n")
  best_res_row = subset(best_res_tbl, V3==pct)
  if(nrow(best_res_row) > 0){
    best_res=as.numeric(subset(best_res_row, V2 == "resolution")$V1)
  }else{
    best_res=0
  }
  
  pct_res_files<-subset(res_files, celltype == pct)
  #cat("  nrow(meta)=", nrow(meta), "\n")

  reductions_rds = file.path(sub_dir, "details", paste0(outFile, ".", celltype_to_filename(pct), ".reductions.rds"))
  reductions<-readRDS(reductions_rds)
  cells<-rownames(reductions$umap)
  subobj<-subset(obj, cells=cells)
  subobj@reductions<-reductions
  
  subumap<-as.data.frame(reductions$subumap@cell.embeddings)
  subumap<-subumap[rownames(subumap) %in% colnames(obj),,drop=FALSE]
  obj@reductions$subumap@cell.embeddings[rownames(subumap), "SUBUMAP_1"] = subumap$subumap_1
  obj@reductions$subumap@cell.embeddings[rownames(subumap), "SUBUMAP_2"] = subumap$subumap_2

  cell_filename = paste0(outFile, ".pre.", celltype_to_filename(pct), ".cell.png")
  save_highlight_cell_plot(cell_filename, subobj, group.by="orig.ident", reduction="subumap", reorder = FALSE) 
  rm(subobj)
  
  if((best_res == 0) | (nrow(pct_res_files) == 0) | (!best_res %in% pct_res_files$resolution)){
    #no corresponding files, which means only one sub cluster
    cat("  only one subcluster\n")
    meta[cells, seurat_clusters] = cluster_index
    meta[cells, resolution_col] = 0
    cluster_index = cluster_index + 1
    
    if(nrow(best_res_row) > 0){
      rename_row = best_res_row[best_res_row$V2 == 0,]
      if(nrow(rename_row) > 0){
        meta[cells, cur_layer] = rename_row$V1[1]
      }
    }
    
    if(output_heatmap){
      allmarkers=c(allmarkers, unlist(ct_top10_map[pct]))
    }
    
    next
  }
  
  cur_res_files = subset(pct_res_files, resolution==best_res)
  file_map = split(cur_res_files$file, cur_res_files$type)
  
  #cat(file_map, "\n")
  
  meta_rds = file_map$meta
  if(!file.exists(meta_rds)){
    stop(meta_rds)
  }
  all_meta<-readRDS(meta_rds)
  all_meta$seurat_clusters_str<-as.character(all_meta$seurat_clusters)
  ncluster=length(unique(all_meta$seurat_clusters))
  
  cat("  best resolution", best_res, "with", ncluster, "clusters\n")
  
  cur_meta<-all_meta[rownames(all_meta) %in% colnames(obj),,drop=F]
  cur_meta$seurat_clusters=as.numeric(as.character(cur_meta$seurat_clusters))
  
  ct_tbl=subset(best_res_row, V2 != "resolution")
  if(nrow(ct_tbl) > 0){
    ct_tbl$V1<-as.character(ct_tbl$V1)
    ct_tbl$V2<-as.character(ct_tbl$V2)
    
    valid=(ct_tbl$V1 %in% all_meta$seurat_clusters_str) | (ct_tbl$V2 %in% all_meta$seurat_clusters_str)
    invalid=ct_tbl[!valid,,drop=F]
    
    if(nrow(invalid) > 0){
      print(invalid)
      stop("Either first column be cluster for merge, or second column be cluster for rename. ")
    }
    
    rename_tbl=ct_tbl[ct_tbl$V2 %in% cur_meta$seurat_clusters_str,,drop=F]
    if(nrow(rename_tbl) > 0){
      for(idx in c(1:nrow(rename_tbl))){
        sc=rename_tbl$V2[idx]
        scname=rename_tbl$V1[idx]
        cur_meta[cur_meta$seurat_clusters_str==sc, "cur_layer"] = scname
      }
    }
    
    merge_tbl=ct_tbl[ct_tbl$V1 %in% cur_meta$seurat_clusters_str,,drop=F]
    if(nrow(merge_tbl) > 0){
      max_index = max(cur_meta$seurat_clusters) + 1
      
      mname=unique(merge_tbl$V2)[1]
      for (mname in unique(merge_tbl$V2)){
        mcts=merge_tbl$V1[merge_tbl$V2==mname]
        cur_meta[cur_meta$seurat_clusters_str %in% mcts, "cur_layer"] = mname
        if (mname == "DELETE"){
          cur_meta[cur_meta$seurat_clusters_str %in% mcts, "seurat_clusters"] = -10000
        }else{
          cur_meta[cur_meta$seurat_clusters_str %in% mcts, "seurat_clusters"] = max_index
          max_index = max_index + 1
        }
      }
      
      valid_clusters=cur_meta$seurat_clusters >= 0
      ctbl<-table(cur_meta$seurat_clusters[valid_clusters])
      ctbl<-ctbl[order(ctbl, decreasing = T)]
      newnames=c(0:(length(ctbl)-1))
      names(newnames)<-names(ctbl)
      cur_meta$seurat_clusters[valid_clusters] = newnames[as.character(cur_meta$seurat_clusters[valid_clusters])]
    }
  }
  
  cur_meta$sub_clusters<-cur_meta$seurat_clusters
  cur_meta$seurat_clusters<-cur_meta$seurat_clusters + cluster_index
  cluster_index = max(cur_meta$seurat_clusters) + 1
  
  meta[rownames(cur_meta), resolution_col] = best_res
  meta[rownames(cur_meta), seurat_clusters] = cur_meta$seurat_clusters
  meta[rownames(cur_meta), cur_layer] = cur_meta$cur_layer
  meta[rownames(cur_meta), "sub_clusters"] = cur_meta$sub_clusters
  
  if(output_heatmap){
    markers_file = file_map$markers
    cur_markers=read.csv(markers_file, header=T, row.names=1)
    cur_top10 = get_top10_markers(cur_markers)
    allmarkers=c(allmarkers, cur_top19$gene)
  }
}

setwd(cur_folder)

obj@meta.data<-meta

if(any(obj$seurat_clusters<0)){
  #there are cells deleted
  cells<-colnames(obj)[obj$seurat_clusters>=0]
  obj<-subset(obj, cells=cells)
}

meta = obj@meta.data
meta[,seurat_cur_layer] = paste0(meta$seurat_clusters, ": ", meta[,cur_layer])
ct<-meta[!duplicated(meta$seurat_cluster),]
ct<-ct[order(ct$seurat_cluster),]

meta[,cur_layer] =factor(meta[,cur_layer], levels=unique(ct[,cur_layer]))
meta[,seurat_cur_layer] =factor(meta[,seurat_cur_layer], levels=ct[,seurat_cur_layer])

obj@meta.data<-meta

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

if(output_heatmap){
  allmarkers<-unique(allmarkers)
  obj<-myScaleData(obj, allmarkers, "RNA")
}

cat("saving final object ...\n")
saveRDS(obj, paste0(outFile, ".final.rds"))

cat("output figures ...\n")

setwd(tmp_folder)

obj$display_layer<-paste0(obj$sub_clusters, ": ", unlist(obj[[cur_layer]]))
gdot=get_bubble_plot(obj, NULL, NULL, bubblemap_file, assay="RNA", orderby_cluster=TRUE, rotate.title=TRUE, group.by=seurat_cur_layer)
galldata<-gdot$data

pre_layer<-unlist(obj[[previous_layer]])

pcts<-unique(pre_layer)
pct<-pcts[1]
for(pct in pcts){
  cells<-colnames(obj)[pre_layer == pct]
  subobj<-subset(obj, cells=cells)
  
  g<-get_dim_plot(obj = subobj, 
                  group.by = "sub_clusters", 
                  label.by = "display_layer", 
                  label=T, 
                  title=pct, 
                  legend.title="", 
                  reduction="subumap", 
                  ncol=1, 
                  random_colors=TRUE)
  
  g<-g+theme(text = element_text(size = 20))
  png(paste0(outFile,".", celltype_to_filename(pct), ".umap.png"), width=2500, height=2000, res=300)
  print(g)
  dev.off()
  
  g_data<-galldata[galldata$id %in% unlist(subobj[[seurat_cur_layer]]),]
  gdot$data<-g_data
  png(paste0(outFile, ".", celltype_to_filename(pct), ".dot.png"), width=bubble_width, height=get_dot_height_vec(g_data$id),res=300)
  print(gdot)
  dev.off()
}

output_celltype_figures(
  obj = obj,
  cell_identity = cur_layer,
  prefix = prefix,
  bubblemap_file = bubblemap_file,
  cell_activity_database = NULL,
  combined_ct_source = NULL,
  group.by="orig.ident",
  name="sample",
  umap_width=2600,
  dot_width = bubble_width,
  cell_identity_order=NULL,
  all_umap = "umap",
  cell_identity_umap = "subumap")

output_celltype_figures(
  obj = obj, 
  cell_identity = seurat_cur_layer, 
  prefix = prefix, 
  bubblemap_file = bubblemap_file,
  cell_activity_database = NULL,
  combined_ct_source = NULL, 
  group.by="orig.ident", 
  name="sample",
  umap_width=3200,
  dot_width = bubble_width,
  cell_identity_order="seurat_clusters",
  all_umap = "umap",
  cell_identity_umap = "subumap")

setwd(cur_folder)

write.csv(obj[["umap"]]@cell.embeddings, paste0(outFile, ".umap.csv"))
write.csv(obj[["subumap"]]@cell.embeddings, paste0(outFile, ".subumap.csv"))

nclusters<-length(unique(obj$seurat_clusters))

if(output_heatmap){
  width<-max(3000, min(10000, nclusters * 150 + 1000))
  height<-max(3000, min(20000, length(allmarkers) * 60 + 1000))
  
  g<-MyDoHeatMap(obj, max_cell=5000, assay="RNA", features = allmarkers, group.by = seurat_cur_layer, angle = 90) + NoLegend()
  png(paste0(prefix, ".top10.heatmap.png"), width=width, height=height, res=300)
  print(g)
  dev.off()
}

g<-get_dim_plot(obj, group.by = "seurat_clusters", label.by=seurat_cur_layer, label.size = 8, legend.title="") + 
  theme(legend.text = element_text(size = 20)) + ggtitle("") 
png(paste0(prefix, ".umap.png"), width=4400, height=2000, res=300)
print(g)
dev.off()

if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
  g<-get_bubble_plot(obj, "seurat_clusters", cur_layer, bubblemap_file, assay="RNA", orderby_cluster=TRUE)
  png(paste0(prefix, ".dot.png"), width=bubble_width, height=get_dot_height(obj, "seurat_clusters"), res=300)
  print(g)
  dev.off()
}

write.csv(table(obj$cell_type, obj$orig.ident), paste0(outFile, ".ct_orig.ident.csv"))
if(!all(obj$orig.ident == obj$sample)){
  write.csv(table(obj$cell_type, obj$sample), paste0(outFile, ".ct_sample.csv"))
}

#obj<-readRDS("mouse_8363.final.rds")

cat("FindAllMarkers for cell type ...\n")

Idents(obj)<-"cell_type"
markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
markers=markers[markers$p_val_adj < 0.05,]
write.csv(markers, paste0(outFile, ".markers.csv"))
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])

top10genes=unique(top10$gene)
obj<-myScaleData(obj, top10genes, assay="RNA")
height<-max(3000, min(20000, length(top10genes) * 60 + 1000))

nclusters<-length(unique(obj$cell_type))
width<-max(3000, min(10000, nclusters * 300 + 1000))
g<-MyDoHeatMap(obj, max_cell = 5000, assay="RNA", features = top10genes, group.by = "cell_type", angle = 90) + NoLegend()
png(paste0(prefix, ".cell_type.top10.heatmap.png"), width=width, height=height, res=300)
print(g)
dev.off()

nclusters<-length(unique(obj$seurat_cell_type))
width<-max(3000, min(10000, nclusters * 300 + 1000))
g<-MyDoHeatMap(obj, max_cell = 5000, assay="RNA", features = top10genes, group.by = "seurat_cell_type", angle = 90) + NoLegend()
png(paste0(prefix, ".seurat_cell_type.top10.heatmap.png"), width=width, height=height, res=300)
print(g)
dev.off()

cat("done ...\n")
