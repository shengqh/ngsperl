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
seurat_clusters = "seurat_clusters"
seurat_cur_layer=paste0("seurat_", cur_layer)
resolution_col = "resolution"

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
  obj<-factorize_layer(obj, previous_layer)
  Idents(obj)<-previous_layer
}

obj<-AddMetaData(obj, obj[[previous_layer]], col.name = cur_layer)
obj<-unfactorize_layer(obj, cur_layer)
obj<-AddMetaData(obj, -1, col.name = seurat_clusters)
obj<-AddMetaData(obj, -1, col.name = resolution_col)

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

best_res_tbl<-read.table(parSampleFile3, sep="\t", header=F)
best_res_tbl$V1 = as.character(best_res_tbl$V1)
best_res_tbl$V2 = as.character(best_res_tbl$V2)

res_files<-read.csv(parFile4, header=T)
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
    die(x)
  }))
}

if(!all(best_res_tbl$V3 %in% res_files$celltype)){
  defined_ct<-unique(best_res_tbl$V3[!(best_res_tbl$V3 %in% res_files$celltype)])
  stop(paste0("those cell types were defined at celltype_subclusters but not exists:", paste0(defined_ct, sep=",")))
}

if(!all(res_files$celltype %in% best_res_tbl$V3)){
  miss_ct<-unique(res_files$celltype[!(res_files$celltype %in% best_res_tbl$V3)])
  stop(paste0("those cell types were not defined at celltype_subclusters :", paste0(miss_ct, collapse=",")))
}

#remove cell type first
remove_cts<-best_res_tbl$V3[best_res_tbl$V2=="resolution" & best_res_tbl$V1 == "0"]
if(length(remove_cts) > 0){
  obj<-subset(obj, !(unlist(obj[[cur_layer]]) %in% remove_cts))
  obj<-factorize_layer(obj, cur_layer)
}

#find markers for cell types
ct_markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
ct_top10<-get_top10_markers(ct_markers)
ct_top10_map<-split(ct_top10$gene, ct_top10$cluster)

layer4map<-split(tiers$Layer4, tiers$Celltype.name)

meta = obj@meta.data
curprefix = prefix

celltypes<-unlist(meta[[previous_layer]])
tblct<-table(celltypes)
tblct<-tblct[order(tblct, decreasing = T)]

previous_celltypes<-names(tblct)
#previous_celltypes<-c("B cells")

DefaultAssay(obj)<-assay

allmarkers<-NULL
cluster_index=0
pct<-previous_celltypes[5]
for(pct in previous_celltypes){
  cat(pct, "\n")
  cells<-rownames(meta)[meta[,previous_layer] == pct]
  best_res_row = subset(best_res_tbl, V3==pct)
  if(nrow(best_res_row) > 0){
    best_res=as.numeric(subset(best_res_row, V2 == "resolution")$V1)
  }else{
    best_res=0.01
  }
  
  if(best_res == 0){
    next
  }
  
  pct_res_files<-subset(res_files, celltype == pct)
  if((nrow(pct_res_files) == 0) | (!best_res %in% pct_res_files$resolution)){
    #no corresponding files, which means only one sub cluster
    cat("  only one subcluster\n")
    meta[cells, seurat_clusters] = cluster_index
    meta[cells, resolution_col] = 0.01
    cluster_index = cluster_index + 1

    allmarkers=c(allmarkers, unlist(ct_top10_map[pct]))
    next
  }
  
  cur_res_files = subset(pct_res_files, resolution==best_res)
  file_map = unlist(split(cur_res_files$file, cur_res_files$type))

  #cat(file_map, "\n")

  meta_rds = file_map['meta']
  if(!file.exists(meta_rds)){
    stop(meta_rds)
  }
  cur_meta<-readRDS(meta_rds)
  cur_meta$seurat_clusters=as.numeric(as.character(cur_meta$seurat_clusters))

  ncluster=length(unique(cur_meta$seurat_clusters))
  cat("  best resolution", best_res, "with", ncluster, "clusters\n")

  ct_tbl=subset(best_res_row, V2 != "resolution")
  if(nrow(ct_tbl) > 0){
    cur_meta$seurat_clusters_str<-as.character(cur_meta$seurat_clusters)
    ct_tbl$V1<-as.character(ct_tbl$V1)
    ct_tbl$V2<-as.character(ct_tbl$V2)
    
    valid=(ct_tbl$V1 %in% cur_meta$seurat_clusters_str) | (ct_tbl$V2 %in% cur_meta$seurat_clusters_str)
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
      for (mname in unique(merge_tbl$V2)){
        mcts=merge_tbl$V1[merge_tbl$V2==mname]
        cur_meta[cur_meta$seurat_clusters_str %in% mcts, "cur_layer"] = mname
        cur_meta[cur_meta$seurat_clusters_str %in% mcts, "seurat_clusters"] = max_index
        max_index = max_index + 1
      }
      
      ctbl<-table(cur_meta$seurat_clusters)
      ctbl<-ctbl[order(ctbl, decreasing = T)]
      newnames=c(0:(length(ctbl)-1))
      names(newnames)<-names(ctbl)
      cur_meta$seurat_clusters = newnames[as.character(cur_meta$seurat_clusters)]
    }
  }

  cur_meta$seurat_clusters<-cur_meta$seurat_clusters + cluster_index
  cluster_index = max(cur_meta$seurat_clusters) + 1

  meta[rownames(cur_meta), resolution_col] = best_res
  meta[rownames(cur_meta), seurat_clusters] = cur_meta$seurat_clusters
  meta[rownames(cur_meta), cur_layer] = cur_meta$cur_layer

  markers_file = file_map['markers']  
  cur_markers=read.csv(markers_file, header=T, row.names=1)
  cur_top19 = get_top10_markers(cur_markers)
  allmarkers=c(allmarkers, cur_top19$gene)
}

meta[,seurat_cur_layer] = paste0(meta$seurat_clusters, ": ", meta[,cur_layer])
ct<-meta[!duplicated(meta$seurat_cluster),]
ct<-ct[order(ct$seurat_cluster),]

meta[,cur_layer] =factor(meta[,cur_layer], levels=unique(ct[,cur_layer]))
meta[,seurat_cur_layer] =factor(meta[,seurat_cur_layer], levels=ct[,seurat_cur_layer])

obj@meta.data<-meta

obj<-subset(obj, seurat_clusters != -1)

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

allmarkers<-unique(allmarkers)
obj<-myScaleData(obj, allmarkers, "RNA")
saveRDS(obj, paste0(outFile, ".obj.rds"))

nclusters<-length(unique(obj$seurat_clusters))

width<-max(3000, min(10000, nclusters * 150 + 1000))
height<-max(3000, min(20000, length(allmarkers) * 60 + 1000))

g<-DoHeatmap(obj, assay="RNA", features = allmarkers, group.by = seurat_cur_layer, angle = 90) + NoLegend()
png(paste0(prefix, ".top10.heatmap.png"), width=width, height=height, res=300)
print(g)
dev.off()

get_cluster_colors<-function(n){
  return(hue_pal()(n))
}

ccolors<-get_cluster_colors(nclusters)

set.seed(random.seed)
scolors<-sample(ccolors, size=nclusters)
g<-DimPlot(obj, group.by = "seurat_clusters", label=T) + ggtitle(cur_layer)+
      scale_color_manual(values=scolors, labels = ct[,seurat_cur_layer])
if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
  g<-g+get_bubble_plot(obj, "seurat_clusters", cur_layer, bubblemap_file, assay="RNA", TRUE)
  g<-g+plot_layout(ncol = 2, widths = c(4, 6))
  width=11000
}else{
  width=4300
}

g<-g+theme(text = element_text(size = 20)) 
png(paste0(prefix, ".umap.png"), width=width, height=4000, res=300)
print(g)
dev.off()
