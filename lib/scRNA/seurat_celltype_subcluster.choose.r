rm(list=ls()) 
outFile='combined'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3='fileList3.txt'
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn/result/combined.final.rds'
parFile2='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_2_subcluster/result/combined.meta.rds'
parFile3='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/essential_genes/result/combined.txt'
parFile4='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_2_subcluster/result/combined.files.csv'

setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result')

### Parameter setting end ###

source("scRNA_func.r")
library(tools)
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
pca_dims=as.numeric(myoptions$pca_dims)

previous_layer<-myoptions$celltype_layer
output_layer<-myoptions$output_layer
seurat_clusters = "seurat_clusters"
seurat_output_layer=paste0("seurat_", output_layer)
resolution_col = "resolution"

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1
essential_genes<-""

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = myoptions$species,
                             curated_markers_file = myoptions$curated_markers_file,
                             HLA_panglao5_file = myoptions$HLA_panglao5_file,
                             layer="Layer4",
                             remove_subtype_str = "",
                             combined_celltype_file = NULL)

cell_activity_database<-ctdef$cell_activity_database


bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)
bubble_width <-ifelse(is.null(myoptions$bubble_width), 6000, as.numeric(myoptions$bubble_width))

prefix<-outFile

if(!exists("obj")){
  obj<-read_object(parFile1, parFile2)

  renamed_layer=paste0(previous_layer, "_renamed")
  if(renamed_layer %in% colnames(obj@meta.data)){
    previous_layer = renamed_layer
  }

  obj<-factorize_layer(obj, previous_layer)
  Idents(obj)<-previous_layer
}
# obj<-readRDS("crs.final.rds")

pre_choose_file=paste0(prefix, ".pre.umap.png")

post_rename_umap = file.path(dirname(parFile4), paste0(outFile, ".post_rename.umap.png"))
if(file.exists(post_rename_umap)){
  file.copy(post_rename_umap, pre_choose_file, overwrite=TRUE)
}else{
  obj<-build_dummy_cluster(obj, label.by=previous_layer, new_cluster_name = "old_clusters")
  g<-get_dim_plot(obj, group.by ="old_clusters", label.by="old_clusters_label", label.size = 8, legend.title="") + 
    theme(legend.text = element_text(size = 20)) + ggtitle("")
  ggsave(pre_choose_file, g, width=3200, height=2000, units="px", dpi=300, bg="white")
  obj$old_clusters<-NULL
  obj$old_clusters_label<-NULL
}

obj<-AddMetaData(obj, obj[[previous_layer]], col.name = output_layer)
obj<-unfactorize_layer(obj, output_layer)
obj<-AddMetaData(obj, -1, col.name = seurat_clusters)
obj<-AddMetaData(obj, -1, col.name = resolution_col)
obj<-AddMetaData(obj, 0, col.name = "sub_clusters")

#set the subuamp using default umap
#it would be replaced cell type by cell type using real sub umap in subclustering
subumap<-Embeddings(obj, reduction = "umap")
colnames(subumap)<-c("SUBUMAP_1", "SUBUMAP_2")
obj[["subumap"]] <- CreateDimReducObject(embeddings = subumap, key = "SUBUMAP_", assay = DefaultAssay(obj))

clcounts<-table(obj[[output_layer]])

if(has_bubblemap){
  allgenes<-rownames(obj)
  genes_df <- read_bubble_genes(bubblemap_file, allgenes, species = myoptions$species)
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
  cells<-colnames(obj)[!(unlist(obj[[output_layer]]) %in% remove_cts)]
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
meta[, output_layer] = as.character(meta[, output_layer])

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

#meta = obj@meta.data
#meta[, cur_layer] = as.character(meta[, cur_layer])

process_keep=function(ct_tbl, cur_meta){
  keep_tbl=ct_tbl |> dplyr::filter(V2 == "KEEP")
  if(nrow(keep_tbl) > 0){
    keep_formula=keep_tbl$V1[1]
    for(keep_formula in keep_tbl$V1) {
      keep_parts=unlist(strsplit(keep_formula, ":"))
      
      keep_cluster=keep_parts[1]
      keep_column=keep_parts[2]
      keep_celltypes_str=keep_parts[3]
      keep_celltypes=unlist(strsplit(keep_celltypes_str, ","))

      cat("  keep cluster:", keep_cluster, ", column:", keep_column, ", celltypes:", keep_celltypes_str, "\n")
      if(!keep_column %in% colnames(cur_meta)){
        stop(paste0("column ", keep_column, " not exists"))
      }

      is_delete = cur_meta$seurat_clusters_str==keep_cluster & (!cur_meta[,keep_column] %in% keep_celltypes)
      delete_meta=cur_meta[is_delete,,drop=F]
      delete_cts=unique(delete_meta[,keep_column])

      delete_cells=sum(is_delete)
      cur_meta$seurat_clusters[is_delete] = -10000
      cat("     deleted cells:", delete_cells, "of", paste0(delete_cts, collpase=","), "\n")
    }
  }
  return(cur_meta)
}

process_naming=function(name_tbl, cur_meta){
  if(nrow(name_tbl) > 0){
    for(idx in c(1:nrow(name_tbl))){
      sc=name_tbl$V2[idx]
      scname=name_tbl$V1[idx]
      cat("  name cluster:", sc, "to", scname, "\n")
      if (scname == "DELETE") {
        is_delete = cur_meta$seurat_clusters_str==sc
        delete_cells=sum(is_delete)
        cur_meta$seurat_clusters[is_delete] = -10000
        cat("     deleted cells:", delete_cells, "of", paste0(delete_cts, collpase=","), "\n")
      }else{
        cur_meta[cur_meta$seurat_clusters_str==sc, "cur_layer"] = scname
      }
    }
  }
  return(cur_meta)
}

meta$seurat_clusters=-1
cluster_index=0
pct<-previous_celltypes[1]
for(pct in previous_celltypes){
  cat(pct, "\n")

  cells<-rownames(meta)[meta[,previous_layer] == pct]

  best_res_row = subset(best_res_tbl, V3==pct)
  if(nrow(best_res_row) > 0){
    best_res_str=subset(best_res_row, V2 == "resolution")$V1
    best_res=as.numeric(best_res_str)
    if(is.na(best_res)){
      best_res=best_res_str
    }
  }else{
    best_res=0
  }
  
  pct_res_files<-subset(res_files, celltype == pct)
  #cat("  nrow(meta)=", nrow(meta), "\n")

  reductions_rds = file.path(sub_dir, "details", paste0(outFile, ".", celltype_to_filename(pct), ".reductions.rds"))
  reductions<-readRDS(reductions_rds)
  
  subumap<-as.data.frame(reductions$subumap@cell.embeddings)
  subumap<-subumap[rownames(subumap) %in% colnames(obj),,drop=FALSE]
  obj@reductions$subumap@cell.embeddings[rownames(subumap), "SUBUMAP_1"] = subumap$subumap_1
  obj@reductions$subumap@cell.embeddings[rownames(subumap), "SUBUMAP_2"] = subumap$subumap_2

  # subobj<-subset(obj, cells=cells)
  # subobj@reductions<-reductions
  # cell_filename = paste0(outFile, ".pre.", celltype_to_filename(pct), ".cell.png")
  # save_highlight_cell_plot(cell_filename, subobj, group.by="orig.ident", reduction="subumap", reorder = FALSE) 
  # rm(subobj)
  
  is_single_cluster=FALSE
  if((best_res == 0) | (nrow(pct_res_files) == 0)){
    is_single_cluster=TRUE
  }
  if(is.numeric(best_res) & (!best_res %in% pct_res_files$resolution)){
    is_single_cluster=TRUE
  }

  if(is_single_cluster){
    #no corresponding files, which means only one sub cluster
    cat("  only one subcluster\n")
    meta[cells, "sub_clusters"] = 0
    meta[cells, seurat_clusters] = cluster_index
    meta[cells, resolution_col] = "0"
    cluster_index = cluster_index + 1
    
    if(nrow(best_res_row) > 0){
      rename_row = best_res_row[best_res_row$V2 == 0,]

      ct_tbl=subset(best_res_row, V2 != "resolution")
      cur_meta=meta[cells,]
      cur_meta$seurat_clusters_str="0"

      if(nrow(rename_row) > 0){
        meta[cells, output_layer] = rename_row$V1[1]
      }

      cur_meta=process_keep(ct_tbl, cur_meta)

      meta[rownames(cur_meta), "seurat_clusters"]=cur_meta$seurat_clusters
    }
    
    if(output_heatmap){
      allmarkers=c(allmarkers, unlist(ct_top10_map[pct]))
    }
    
    print(table(meta[,output_layer], meta$seurat_clusters))

    next
  }

  if(is.character(best_res)){
    #assign cell type based on SignacX/SingleR/Azimuth
    cat("  assign cluster by", best_res, "\n")
    stopifnot(all(cells %in% rownames(meta)))

    cur_meta=meta[cells,]
    cur_meta$cur_layer = cur_meta[,best_res]
    cur_meta$cur_layer[is.na(cur_meta$cur_layer)] = "DELETE"

    ct_tbl=subset(best_res_row, V2 != "resolution")
    if(!"DELETE" %in% ct_tbl$V2){
      #The cell types not in list will be deleted.
      ct_tbl=rbind(ct_tbl, data.frame(V1="OTHERS", V2="DELETE", V3=pct))
    }

    source_cts=unique(ct_tbl$V1)
    anno_cts=setdiff(source_cts, "OTHERS")

    source_ct=source_cts[2]
    for(source_ct in source_cts){
      target_ct=ct_tbl$V2[ct_tbl$V1==source_ct]
      cat("   ", source_ct, "to", target_ct, "\n")
      if(source_ct == "OTHERS"){
        cur_cells=rownames(cur_meta)[!cur_meta[,best_res] %in% anno_cts]
      }else{
        cur_cells=rownames(cur_meta)[cur_meta[,best_res] %in% source_ct]
      }
      cur_meta[cur_cells, "cur_layer"] = target_ct
    }
    cur_meta$seurat_clusters[cur_meta$cur_layer == "DELETE"] = -10000

    tbl=table(cur_meta$cur_layer)
    tbl=tbl[names(tbl) != "DELETE"]
    tbl=tbl[order(tbl, decreasing = T)]

    cur_index = 0
    for(tname in names(tbl)){
      cur_meta$seurat_clusters[cur_meta$cur_layer == tname] = cur_index
      cur_index = cur_index + 1
    }
  
    cur_meta$sub_clusters<-cur_meta$seurat_clusters
    cur_meta$seurat_clusters<-cur_meta$seurat_clusters + cluster_index
    cluster_index = max(cur_meta$seurat_clusters) + 1
    
    meta[rownames(cur_meta), resolution_col] = best_res
    meta[rownames(cur_meta), seurat_clusters] = cur_meta$seurat_clusters
    meta[rownames(cur_meta), output_layer] = cur_meta$cur_layer
    meta[rownames(cur_meta), "sub_clusters"] = cur_meta$sub_clusters

    if(output_heatmap){
      allmarkers=c(allmarkers, unlist(ct_top10_map[pct]))
    }
    
    print(table(meta[,output_layer], meta$seurat_clusters))

    next
  }

  cur_res_files = subset(pct_res_files, resolution==best_res)
  file_map = split(cur_res_files$file, cur_res_files$type)
  
  #cat(file_map, "\n")
  
  meta_rds = file_map$meta
  if(!file.exists(meta_rds)){
    stop(paste0("missing meta file:", meta_rds))
  }
  all_meta<-readRDS(meta_rds)
  all_meta$seurat_clusters_str<-as.character(all_meta$seurat_clusters)
  ncluster=length(unique(all_meta$seurat_clusters))

  #Assigned cell type is saved in "cur_layer" column of the meta data table
  
  cat("  best resolution", best_res, "with", ncluster, "clusters\n")
  
  cur_meta<-all_meta[rownames(all_meta) %in% colnames(obj),,drop=F]
  cur_meta$seurat_clusters=as.numeric(as.character(cur_meta$seurat_clusters))
  
  ct_tbl=subset(best_res_row, V2 != "resolution")
  if(nrow(ct_tbl) > 0){
    ct_tbl$V1<-as.character(ct_tbl$V1)
    ct_tbl$V2<-as.character(ct_tbl$V2)

    move_tbl=ct_tbl |> dplyr::filter(V2 == "MOVE")
    if(nrow(move_tbl) > 0){
      move_formula=move_tbl$V1[1]
      for(move_formula in move_tbl$V1) {
        move_parts=unlist(strsplit(move_formula, ":"))
        
        move_cluster=move_parts[1]
        cur_name=unique(cur_meta$cur_layer[cur_meta$seurat_clusters_str==move_cluster])
        move_column=move_parts[2]
        move_celltypes_str=move_parts[3]
        move_celltypes=unlist(strsplit(move_celltypes_str, ","))
        to_cluster=as.numeric(move_parts[4])

        cat("  moving", move_celltypes_str, "annotated by", move_column, "from cluster", move_cluster, "to", to_cluster, "\n")
        if(!move_column %in% colnames(cur_meta)){
          stop(paste0("column ", move_column, " not exists"))
        }

        if(move_cluster == -1){
          is_move = cur_meta[,move_column] %in% move_celltypes
        }else{
          is_move = cur_meta$seurat_clusters_str==move_cluster & cur_meta[,move_column] %in% move_celltypes
        }

        move_cells=sum(is_move)
        cur_meta$seurat_clusters[is_move] = to_cluster
        cat("     moved cells:", move_cells, "\n")

        cur_meta$seurat_clusters_str<-as.character(cur_meta$seurat_clusters)        
      }
    }

    name_or_merge_tbl=ct_tbl |> 
      dplyr::filter(!(V2 %in% c("DELETE", "KEEP", "RENAME", "MOVE")))
    
    if(nrow(name_or_merge_tbl) > 0){
      valid=name_or_merge_tbl$V2 %in% cur_meta$seurat_clusters_str | name_or_merge_tbl$V1 %in% cur_meta$seurat_clusters_str
      invalid=name_or_merge_tbl[!valid,,drop=F]
      
      if(nrow(invalid) > 0){
        print(invalid)
        stop("Second column should be cluster id for naming and first column should be cluster id for merging. ")
      }
    }

    name_tbl=name_or_merge_tbl |> dplyr::filter(V2 %in% cur_meta$seurat_clusters_str)
    cur_meta = process_naming(name_tbl, cur_meta)

    delete_by_cluster_tbl=ct_tbl |> dplyr::filter(V2 == "DELETE") |> dplyr::filter(V1 %in% cur_meta$seurat_clusters_str)
    if(nrow(delete_by_cluster_tbl) > 0){
      delete_scs = delete_by_cluster_tbl$V1
      cat("  deleting cluster", paste0(delete_scs, collapse=","), "of", pct, "\n")

      id_delete = cur_meta$seurat_clusters_str %in% delete_scs
      delete_cells=sum(id_delete)
      cur_meta$seurat_clusters[id_delete] = -10000
      cat("     deleted", delete_cells, "outof", nrow(cur_meta), "cells\n")
    }

    delete_by_formula_tbl=ct_tbl |> dplyr::filter(V2 == "DELETE") |> dplyr::filter(!(V1 %in% cur_meta$seurat_clusters_str))
    if(nrow(delete_by_formula_tbl) > 0){
      delete_formula=delete_by_formula_tbl$V1[1]
      for(delete_formula in delete_by_formula_tbl$V1) {
        delete_parts=unlist(strsplit(delete_formula, ":"))
        
        delete_cluster=delete_parts[1]
        cur_name=unique(cur_meta$cur_layer[cur_meta$seurat_clusters_str==delete_cluster])
        delete_column=delete_parts[2]
        delete_celltypes_str=delete_parts[3]
        delete_celltypes=unlist(strsplit(delete_celltypes_str, ","))

        cat("  deleting", delete_celltypes_str, "annotated by", delete_column, "from cluster", delete_cluster, "(", cur_name, ")\n")
        if(!delete_column %in% colnames(cur_meta)){
          stop(paste0("column ", delete_column, " not exists"))
        }

        is_delete = cur_meta$seurat_clusters_str==delete_cluster & cur_meta[,delete_column] %in% delete_celltypes
        delete_cells=sum(is_delete)
        cur_meta$seurat_clusters[is_delete] = -10000
        cat("     deleted cells:", delete_cells, "\n")
      }
    }

    keep_tbl=ct_tbl |> dplyr::filter(V2 == "KEEP")
    if(nrow(keep_tbl) > 0){
      keep_formula=keep_tbl$V1[1]
      for(keep_formula in keep_tbl$V1) {
        keep_parts=unlist(strsplit(keep_formula, ":"))
        
        keep_cluster=keep_parts[1]
        keep_column=keep_parts[2]
        keep_celltypes_str=keep_parts[3]
        keep_celltypes=unlist(strsplit(keep_celltypes_str, ","))

        cat("  keep cluster:", keep_cluster, ", column:", keep_column, ", celltypes:", keep_celltypes_str, "\n")
        if(!keep_column %in% colnames(cur_meta)){
          stop(paste0("column ", keep_column, " not exists"))
        }

        is_delete = cur_meta$seurat_clusters_str==keep_cluster & (!cur_meta[,keep_column] %in% keep_celltypes)
        delete_meta=cur_meta[is_delete,,drop=F]
        delete_cts=unique(delete_meta[,keep_column])

        delete_cells=sum(is_delete)
        cur_meta$seurat_clusters[is_delete] = -10000
        cat("     deleted cells:", delete_cells, "of", paste0(delete_cts, collpase=","), "\n")
      }
    }
    
    rename_tbl=ct_tbl |> dplyr::filter(V2 == "RENAME")
    if(nrow(rename_tbl) > 0){
      rename_formula=rename_tbl$V1[1]
      for(rename_formula in rename_tbl$V1) {
        rename_parts=unlist(strsplit(rename_formula, ":"))

        if(length(rename_parts) != 4){
          stop(paste0("  rename formula should be cluster:column:source_celltypes:target_celltype, now we get", rename_formula))
        }
        
        rename_cluster=rename_parts[1]
        rename_column=rename_parts[2]
        rename_celltypes_str=rename_parts[3]
        rename_target_celltype=rename_parts[4]
        rename_celltypes=unlist(strsplit(rename_celltypes_str, ","))

        cat("  rename cluster:", rename_cluster, ", column:", rename_column, ", celltypes:", rename_celltypes_str, "to", rename_target_celltype, "\n")
        if(!rename_column %in% colnames(cur_meta)){
          stop(paste0("column ", rename_column, " not exists"))
        }

        if(rename_cluster == -1){ #apply to all clusters
          is_rename = cur_meta[,rename_column] %in% rename_celltypes
        }else{
          is_rename = cur_meta$seurat_clusters_str==rename_cluster & (cur_meta[,rename_column] %in% rename_celltypes)
        }
        rename_cells=sum(is_rename)
        max_seurat_clusters = max(cur_meta$seurat_clusters)
        cur_meta[is_rename, "cur_layer"] = rename_target_celltype
        cur_meta[is_rename, "seurat_clusters"] = max_seurat_clusters + 1

        cat("     renamed cells:", rename_cells, "\n")
      }
    }

    merge_tbl=name_or_merge_tbl |> dplyr::filter(V1 %in% cur_meta$seurat_clusters_str)
    if(nrow(merge_tbl) > 0){
      max_index = max(cur_meta$seurat_clusters) + 1
      
      mname=unique(merge_tbl$V2)[1]
      for (mname in unique(merge_tbl$V2)){
        mcts=merge_tbl$V1[merge_tbl$V2==mname]
        n_cells = sum(cur_meta$seurat_clusters %in% mcts)
        cur_meta[cur_meta$seurat_clusters_str %in% mcts, "cur_layer"] = mname
        if (mname == "DELETE"){
          cur_meta[cur_meta$seurat_clusters_str %in% mcts, "seurat_clusters"] = -10000
          cat("     deleted", n_cells, "cells in", paste0(mcts, collpase=","), "\n")
        }else{
          cur_meta[cur_meta$seurat_clusters_str %in% mcts, "seurat_clusters"] = max_index
          cat("     merged", n_cells, "cells to", mname, "in", paste0(mcts, collpase=","), "\n")
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
  meta[rownames(cur_meta), output_layer] = cur_meta$cur_layer
  meta[rownames(cur_meta), "sub_clusters"] = cur_meta$sub_clusters
  
  if(output_heatmap){
    markers_file = file_map$markers
    cur_markers=read.csv(markers_file, header=T, row.names=1)
    cur_top10 = get_top10_markers(cur_markers)
    allmarkers=c(allmarkers, cur_top19$gene)
  }

  print(table(meta[,output_layer], meta$seurat_clusters))
}

setwd(cur_folder)

obj@meta.data<-meta

if(any(obj$seurat_clusters<0)){
  #there are cells deleted
  cells<-colnames(obj)[obj$seurat_clusters>=0]
  obj<-subset(obj, cells=cells)
}

# fix the seurat_clusters. Since we might delete some subclusters, the index might have gap.
# meta=readRDS(meta_rds)
meta = obj@meta.data
old_clusters = names(table(meta$seurat_clusters))
new_clusters = 0:(length(old_clusters)-1)
map_df=data.frame(old_clusters, new_clusters)
meta$seurat_clusters = map_df$new_clusters[match(meta$seurat_clusters, map_df$old_clusters)]

meta[,seurat_output_layer] = paste0(meta$seurat_clusters, ": ", meta[,output_layer])
ct<-meta[!duplicated(meta$seurat_cluster),]
ct<-ct[order(ct$seurat_cluster),]

meta[,output_layer] =factor(meta[,output_layer], levels=unique(ct[,output_layer]))
meta[,seurat_output_layer] =factor(meta[,seurat_output_layer], levels=ct[,seurat_output_layer])

obj@meta.data<-meta

meta_csv = paste0(outFile, ".meta.csv")
write.csv(obj@meta.data, meta_csv)
meta_rds = paste0(outFile, ".meta.rds")
saveRDS(obj@meta.data, meta_rds)

if(output_heatmap){
  allmarkers<-unique(allmarkers)
  obj<-myScaleData(obj, allmarkers, "RNA")
}

# Keep the original UMAP, don't redo PCA and UMAP since it doesn't work for integration.
# cat("redo PCA ...\n")
# obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)
# cat("redo UMAP ...\n")
# obj <- RunUMAP(object = obj, dims=c(1:pca_dims), verbose = FALSE)

cat("saving final object ...\n")
final_file = paste0(outFile, ".final.rds")
saveRDS(obj, final_file)
#obj=readRDS(final_file)

cat("calculate md5sum ...\n")
md5 = md5sum(final_file)
writeLines(md5, paste0(final_file, ".md5"))

cat("output figures ...\n")

setwd(tmp_folder)

obj$display_layer<-paste0(obj$sub_clusters, ": ", unlist(obj[[output_layer]]))
gdot=get_bubble_plot(obj, 
  NULL, 
  NULL, 
  bubblemap_file, 
  assay="RNA", 
  orderby_cluster=TRUE, 
  rotate.title=TRUE, 
  group.by=seurat_output_layer,
  species=myoptions$species)
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
  ggsave( paste0(outFile,".", celltype_to_filename(pct), ".umap.png"), 
          g,
          width=2500, 
          height=2000, 
          units="px",
          dpi=300,
          bg="white")
  
  g_data<-galldata[galldata$id %in% unlist(subobj[[seurat_output_layer]]),]
  gdot$data<-g_data
  ggsave( paste0(outFile, ".", celltype_to_filename(pct), ".dot.png"), 
          gdot,
          width=get_dot_width(gdot, min_width=bubble_width), 
          height=get_dot_height_vec(g_data$id),
          units="px",
          dpi=300,
          bg="white")
}

if(length(unique(obj@meta.data[,seurat_output_layer])) > 18){
  width=3500
}else{
  width=2600
}
output_celltype_figures(
  obj = obj, 
  cell_identity = seurat_output_layer, 
  prefix = prefix, 
  bubblemap_file = bubblemap_file,
  cell_activity_database = NULL,
  combined_ct_source = NULL, 
  group.by="orig.ident", 
  name="sample",
  umap_width=width,
  dot_width = bubble_width,
  cell_identity_order="seurat_clusters",
  all_umap = "umap",
  cell_identity_umap = "subumap",
  species=myoptions$species)

#if there are multiple clusters with same cell type name, draw cell type level figures
if(length(unlist(unique(obj@meta.data[,output_layer]))) != length(unlist(unique(obj@meta.data[,seurat_output_layer])))) {
  if(length(unique(obj@meta.data[,output_layer])) > 18){
    width=3500
  }else{
    width=2600
  }
  output_celltype_figures(
    obj = obj,
    cell_identity = output_layer,
    prefix = prefix,
    bubblemap_file = bubblemap_file,
    cell_activity_database = NULL,
    combined_ct_source = NULL,
    group.by="orig.ident",
    name="sample",
    umap_width=width,
    dot_width = bubble_width,
    cell_identity_order=NULL,
    all_umap = "umap",
    cell_identity_umap = "subumap",
    species=myoptions$species)
}

setwd(cur_folder)

write.csv(obj[["umap"]]@cell.embeddings, paste0(outFile, ".umap.csv"))
write.csv(obj[["subumap"]]@cell.embeddings, paste0(outFile, ".subumap.csv"))

nclusters<-length(unique(obj$seurat_clusters))

if(output_heatmap){
  g<-MyDoHeatMap(obj, max_cell=5000, assay="RNA", features = allmarkers, group.by = seurat_output_layer, angle = 90) + NoLegend()
  ggsave( paste0(prefix, ".top10.heatmap.png"), 
          g,
          width=get_heatmap_width(nclusters), 
          height=get_heatmap_height(length(allmarkers)),
          units="px",
          dpi=300,
          bg="white")
}

g<-get_dim_plot(obj, group.by = "seurat_clusters", label.by=seurat_output_layer, label.size = 8, legend.title="") + 
  theme(legend.text = element_text(size = 20)) + ggtitle("") 
if(length(unique(obj$seurat_clusters)) > 18){
  width=6000
}else{
  width=4000
}
ggsave( paste0(prefix, ".umap.png"), 
        g,
        width=width, 
        height=2000, 
        units="px", 
        dpi=300, 
        bg="white")

data_norm=get_seurat_average_expression(obj, seurat_output_layer)
predict_celltype<-ORA_celltype(data_norm,cell_activity_database$cellType,cell_activity_database$weight)
saveRDS(predict_celltype, paste0(outFile, ".cta.rds"))
Plot_predictcelltype_ggplot2( predict_celltype, 
                              filename=paste0(outFile, ".cta.png"),
                              is_validation=TRUE)

if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
  g<-get_bubble_plot(obj, 
    "seurat_clusters", 
    output_layer, 
    bubblemap_file, 
    assay="RNA", 
    orderby_cluster=TRUE,
    species=myoptions$species)
  ggsave( paste0(prefix, ".dot.png"), 
          g,
          width=get_dot_width(g, min_width=bubble_width), 
          height=get_dot_height(obj, "seurat_clusters"), 
          units="px", 
          dpi=300, 
          bg="white")
}

write.csv(table(obj$cell_type, obj$orig.ident), paste0(outFile, ".ct_orig.ident.csv"))
if(!all(obj$orig.ident == obj$sample)){
  write.csv(table(obj$cell_type, obj$sample), paste0(outFile, ".ct_sample.csv"))
}

#obj<-readRDS("combined.final.rds")

cat("FindAllMarkers for cell type ...\n")

Idents(obj)<-"cell_type"
markers=FindAllMarkers(obj, assay="RNA", only.pos=TRUE, min.pct=min.pct, logfc.threshold=logfc.threshold)
markers=markers[markers$p_val_adj < 0.05,]

markers_file=paste0(outFile, ".markers.csv")
write.csv(markers, markers_file)

markers<-read.csv(markers_file, header=T, row.names=1)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = .data[["avg_log2FC"]])
top10genes=unique(top10$gene)
obj<-myScaleData(obj, top10genes, assay="RNA")

g<-MyDoHeatMap(obj, max_cell = 5000, assay="RNA", features = top10genes, group.by = "cell_type", angle = 90) + NoLegend()
ggsave( paste0( prefix, ".cell_type.top10.heatmap.png"), 
        g,
        width=get_heatmap_width(length(unique(obj$cell_type))), 
        height=get_heatmap_height(length(top10genes)), 
        units="px",
        dpi=300,
        bg="white",
        limitsize=FALSE)

g<-MyDoHeatMap(obj, max_cell = 5000, assay="RNA", features = top10genes, group.by = "seurat_cell_type", angle = 90) + NoLegend()
ggsave( paste0( prefix, ".seurat_cell_type.top10.heatmap.png"), 
        g,
        width=get_heatmap_width(length(unique(obj$seurat_cell_type))), 
        height=get_heatmap_height(length(top10genes)), 
        units="px",
        dpi=300,
        bg="white",
        limitsize=FALSE)

cat("done ...\n")

