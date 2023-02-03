rm(list=ls()) 
outFile='P9270'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1='/workspace/shengq2/charles_flynn/20230105_9270_scRNA_dog/seurat_sct_merge/result/P9270.final.rds'
parFile2=''
parFile3='/workspace/shengq2/charles_flynn/20230105_9270_scRNA_dog/essential_genes/result/P9270.txt'


setwd('/workspace/shengq2/charles_flynn/20230105_9270_scRNA_dog/seurat_sct_merge_dr0.5_individual/result')

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
library('rmarkdown')

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
remove_subtype_str<-myoptions$remove_subtype
annotate_tcell<-is_one(myoptions$annotate_tcell)
HLA_panglao5_file<-myoptions$HLA_panglao5_file
tcell_markers_file<-myoptions$tcell_markers_file
assay=ifelse(by_sctransform, "SCT", "RNA")
by_harmony<-reduction=="harmony"
regress_by_percent_mt<-is_one(myoptions$regress_by_percent_mt)
redo_harmony<-is_one(myoptions$redo_harmony, 0)
resolution=as.numeric(myoptions$dynamic_by_one_resolution)
curated_markers_file=myoptions$curated_markers_file
by_individual_sample=is_one(myoptions$by_individual_sample)

if(by_individual_sample){
  by_harmony<-FALSE
}

layer=ifelse(is.null(myoptions$layer), "Layer4", myoptions$layer)

if(regress_by_percent_mt){
  vars.to.regress="percent.mt"
}else{
  vars.to.regress=NULL
}

essential_genes=read.table(parFile3, sep="\t" ,header=F)$V1

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is_file_empty(bubblemap_file)

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

if(!exists('parSampleFile4')){
  parSampleFile4=""
}

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = species,
                             curated_markers_file = curated_markers_file,
                             HLA_panglao5_file = HLA_panglao5_file,
                             layer=layer,
                             remove_subtype_str = remove_subtype_str,
                             combined_celltype_file = parSampleFile4)

cell_activity_database<-ctdef$cell_activity_database

layer2map<-ctdef$celltype_map

combined_ct<-ctdef$combined_celltypes
combined_ct_source<-ctdef$combined_celltype_source

replace_cts=layer2map %in% names(combined_ct)
layer2map[replace_cts] = combined_ct[layer2map[replace_cts]]

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
  genes_df <- read_bubble_genes(bubblemap_file, allgenes)
  bubble_genes<-unique(genes_df$gene)
}

obj[["layer0"]]<-"Unassigned"
obj[["layer0_clusters"]]<-0
obj[["layer0_raw"]]<-"Unassigned"

cbind_celltype<-function(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts){
  if(is.null(cur_layermap)){
    return(cur_cts)
  }
  layer_ids<-unlist(cur_layermap[new_cluster_ids])
  names(layer_ids) <- colnames(data_norm)
  
  oldcluster<-subobj[[cluster]][[1]]
  cur_cts$seurat_clusters=oldcluster
  cur_cts$raw_cell_type<-new_cluster_ids[oldcluster]
  cur_cts$raw_seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$raw_cell_type) 
  cur_cts$cell_type<-layer_ids[oldcluster]
  cur_cts$seurat_cell_type<-paste0(cur_cts$seurat_cluster, ": ", cur_cts$cell_type)

  return(cur_cts)
}

if(0){
  previous_layer<-"layer0"
  cur_layer="layer4"
  cur_layermap=layer2map
  previous_celltypes<-unique(obj@meta.data[[previous_layer]])
  iter=1
}

get_empty_files<-function(){
  files = data.frame("previous_layer"=character(),
                     "cur_layer"=character(),
                     "pct"=character(),
                     "type"=character(),
                     "fname"=character())
  return(files)
}

iterate_celltype<-function(obj, 
                           previous_celltypes, 
                           previous_layer, 
                           cur_layer, 
                           cur_layermap, 
                           npcs, 
                           resolution, 
                           random.seed, 
                           by_sctransform, 
                           by_harmony, 
                           curprefix, 
                           iter, 
                           vars.to.regress,
                           bubblemap_file, 
                           essential_genes){
  meta = obj@meta.data
  
  assay=ifelse(by_sctransform, "SCT", "RNA")

  files = get_empty_files()
  
  all_cur_cts<-NULL
  pct<-previous_celltypes[length(previous_celltypes)]
  
  #previous_celltypes<-c("Platelets")
  for(pct in previous_celltypes){
    pct_str = celltype_to_filename(pct)

    key = paste0("iter", iter, ": ", pct, ":")
    cells<-rownames(meta)[meta[,previous_layer] == pct]
    if(length(cells) == 0){#no cell left for this cell type
      next
    }
    
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
      subobj = sub_cluster(subobj = subobj, 
                            assay =  assay, 
                            by_sctransform = by_sctransform, 
                            by_harmony = by_harmony, 
                            redo_harmony = redo_harmony,
                            curreduction = curreduction, 
                            k_n_neighbors = k_n_neighbors,
                            u_n_neighbors = u_n_neighbors,
                            random.seed = random.seed,
                            resolutions = resolution,
                            cur_npcs = cur_npcs, 
                            cur_pca_dims = cur_pca_dims,
                            vars.to.regress = vars.to.regress,
                            essential_genes = essential_genes
                            )
    }else{
      cat(key, "FindNeighbors\n")
      subobj<-FindNeighbors(object=subobj, reduction=curreduction, k.param=k_n_neighbors, dims=cur_pca_dims, verbose=FALSE)

      cat(key, "FindClusters\n")
      subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolution, verbose=FALSE)
    }
    
    cat(key, "Cell type annotation\n")
    cur_cts<-subobj[[previous_layer]]

    cluster<-"seurat_clusters"
    data_norm=get_seurat_average_expression(subobj, cluster)
    
    predict_celltype<-ORA_celltype(data_norm,cell_activity_database$cellType,cell_activity_database$weight)
    
    cta_rds_file=paste0(curprefix, ".", pct_str, ".cta.rds")
    saveRDS(predict_celltype, cta_rds_file)
    files<-rbind(files, c(previous_layer, cur_layer, pct, "cta_rds", cta_rds_file))

    if(length(predict_celltype$max_cta) > 1){
      cta_png_file=paste0(curprefix, ".", pct_str, ".cta.png")
      Plot_predictcelltype( predict_celltype, 
                            filename=cta_png_file)
      files<-rbind(files, c(previous_layer, cur_layer, pct, "cta_png", cta_png_file))
    }

    new_cluster_ids<-names(predict_celltype$max_cta)
    
    cur_cts<-cbind_celltype(subobj, data_norm, cluster, new_cluster_ids, cur_layermap, cur_cts)
    stopifnot(all(colnames(subobj) == rownames(cur_cts)))
    subobj = AddMetaData(subobj, cur_cts$seurat_clusters, "seurat_clusters")
    subobj = AddMetaData(subobj, cur_cts$cell_type, "cell_type")
    subobj = AddMetaData(subobj, cur_cts$seurat_cell_type, "seurat_cell_type")
    subobj = AddMetaData(subobj, cur_cts$raw_cell_type, "raw_cell_type")
    subobj = AddMetaData(subobj, cur_cts$raw_seurat_cell_type, "raw_seurat_cell_type")
    
    #using RNA assay for visualization
    DefaultAssay(subobj)<-assay

    g0<-DimPlot(obj, label=F, cells.highlight =cells) + ggtitle(pct) + scale_color_discrete(type=c("gray", "red"), labels = c("others", pct))
    g1<-DimPlot(subobj, reduction="oumap", group.by = "cell_type", label=T) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("New cell type in old UMAP")
    g2<-get_dim_plot(subobj, reduction="oumap", group.by="seurat_clusters", label.by="raw_seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat raw cell type in old UMAP")
    g3<-get_dim_plot(subobj, reduction="oumap", group.by="seurat_clusters", label.by="seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat cell type in old UMAP")
    g<-g0+g1+g2+g3+plot_layout(nrow=2)
    umap_celltype_file = paste0(curprefix, ".", pct_str, ".old_umap.png")
    png(umap_celltype_file, width=4600, height=4000, res=300)
    print(g)
    dev.off()

    files<-rbind(files, c(previous_layer, cur_layer, pct, "old_umap", umap_celltype_file))

    if(pct != "Unassigned"){
      g1<-get_dim_plot(subobj, group.by="seurat_clusters", label.by="raw_seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Raw cell type in new UMAP")
      g2<-get_dim_plot(subobj, group.by="seurat_clusters", label.by="seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat cell type in new UMAP")

      g<-g1+g2+plot_layout(nrow=1)
      umap_cluster_file = paste0(curprefix, ".", pct_str, ".new_umap.png")
      png(umap_cluster_file, width=4600, height=2000, res=300)
      print(g)
      dev.off()
      files<-rbind(files, c(previous_layer, cur_layer, pct, "new_umap", umap_cluster_file))
    }

    dot_width=4400
    g<-get_bubble_plot(subobj, cur_res=cluster, "raw_cell_type", bubblemap_file, assay="RNA", orderby_cluster = FALSE)
    dot_file = paste0(curprefix, ".", pct_str, ".dot.png")
    png(dot_file, width=dot_width, height=get_dot_height(subobj, cluster), res=300)
    print(g)
    dev.off()
    files<-rbind(files, c(previous_layer, cur_layer, pct, "dot", dot_file))

    all_cur_cts<-rbind(all_cur_cts, cur_cts)
    
    # if(previous_layer == "layer0"){
    #   obj[['umap']] = subobj[['umap']]
    # }
    rm(subobj)
  }
  colnames(files)<-colnames(get_empty_files())
  return(list("all_cur_cts"=all_cur_cts, "files"=files))
}

layer_cluster_celltype<-function(obj, 
                                 previous_layer, 
                                 cur_layer, 
                                 cur_layermap, 
                                 npcs, 
                                 resolution, 
                                 random.seed, 
                                 by_sctransform, 
                                 by_harmony, 
                                 prefix, 
                                 vars.to.regress,
                                 bubblemap_file,
                                 essential_genes){
  meta<-obj@meta.data
  
  previous_celltypes<-unique(meta[[previous_layer]])

  files = get_empty_files()

  if(length(previous_celltypes) > 0){
    cat("cluster and annotate cell types:", paste0(previous_celltypes, ", "))
    iter = 1
    while(TRUE){
      cat("Iteration ", iter, "\n")
      
      iter_name=paste0("iter", iter)

      previous_celltypes<-previous_celltypes[order(previous_celltypes)]
      
      curprefix = paste0(prefix, ".iter", iter)

      iter_meta_file = paste0(curprefix, ".csv")
      iter_meta_rds = paste0(curprefix, ".rds")

      lst<-iterate_celltype(obj, 
                            previous_celltypes, 
                            previous_layer, 
                            iter_name, 
                            cur_layermap, 
                            npcs, 
                            resolution, 
                            random.seed, 
                            by_sctransform, 
                            by_harmony, 
                            curprefix, 
                            iter, 
                            vars.to.regress,
                            bubblemap_file, 
                            essential_genes)

      all_cur_cts<-lst$all_cur_cts
      cur_files<-lst$files

      files<-rbind(files, cur_files)
      
      stopifnot(all(rownames(all_cur_cts) %in% colnames(obj)))
      
      iter_clusters = paste0(iter_name, "_clusters")
      iter_raw = paste0(iter_name, "_raw")

      obj[[iter_name]] = obj[[previous_layer]]
      obj[[iter_clusters]] = obj[[paste0(previous_layer, "_clusters")]]
      obj[[iter_raw]] = obj[[paste0(previous_layer, "_raw")]]

      obj@meta.data[rownames(all_cur_cts), iter_name]=unlist(all_cur_cts[, "cell_type"])
      obj@meta.data[rownames(all_cur_cts), iter_clusters]=unlist(as.character(all_cur_cts[, "seurat_clusters"]))
      obj@meta.data[rownames(all_cur_cts), iter_raw]=unlist(all_cur_cts[, "raw_cell_type"])
      
      pre_disagree<-all_cur_cts[all_cur_cts[, previous_layer] != all_cur_cts[,'cell_type'],,drop=F]
      if(nrow(pre_disagree) > 0){
        cat("Found unmatched cell type\n")
        print(table(pre_disagree[,previous_layer], pre_disagree$cell_type))
        write.csv(obj@meta.data, iter_meta_file)
        saveRDS(obj@meta.data, iter_meta_rds)
        previous_celltypes = unique(c(unlist(pre_disagree[,previous_layer]), unlist(pre_disagree$cell_type)))
        previous_celltypes<-previous_celltypes[previous_celltypes %in% unlist(unique(obj[[iter_name]]))]
        previous_layer = iter_name
        iter = iter + 1
        next
      }else{
        obj[[cur_layer]] = obj[[iter_name]]
        obj[[paste0(cur_layer, "_clusters")]] = obj[[iter_clusters]]
        obj[[paste0(cur_layer, "_raw")]] = obj[[iter_raw]]
        write.csv(obj@meta.data, iter_meta_file)
        saveRDS(obj@meta.data, iter_meta_rds)
        break
      }
    }
  }
  
  #using RNA assay for visualization
  DefaultAssay(obj)<-"RNA"

  g<-DimPlot(obj, group.by = cur_layer, label=T)

  png(paste0(prefix, ".", cur_layer, ".umap.png"), width=3300, height=3000, res=300)
  print(g)
  dev.off()

  if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
    g2<-get_bubble_plot(obj, NA, cur_layer, bubblemap_file, assay="RNA")
    png(paste0(prefix, ".", cur_layer, ".dot.png"), width=4400, height=2000, res=300)
    print(g2)
    dev.off()
  }

  write.csv(obj@meta.data, file=paste0(prefix, ".", cur_layer, ".meta.csv"))
  saveRDS(obj@meta.data, file=paste0(prefix, ".", cur_layer, ".meta.rds"))

  write.csv(all_cur_cts, file=paste0(prefix, ".", cur_layer, ".details.csv"))
  
  return(list(obj=obj, files=files))
}

if(0){
  previous_layer = "layer0"
  cur_layer = "layer4"
  cur_layermap = layer2map
}

do_analysis<-function(tmp_folder,
                      cur_folder,
                      obj, 
                      layer2map, 
                      npcs, 
                      resolution, 
                      random.seed, 
                      by_sctransform, 
                      by_harmony, 
                      prefix, 
                      vars.to.regress, 
                      bubblemap_file, 
                      essential_genes,
                      by_individual_sample ) {
  setwd(tmp_folder)
  reslist1<-layer_cluster_celltype(obj = obj,
                              previous_layer = "layer0", 
                              cur_layer = "layer4", 
                              cur_layermap = layer2map, 
                              npcs = npcs, 
                              resolution = resolution, 
                              random.seed = random.seed, 
                              by_sctransform = by_sctransform, 
                              by_harmony = by_harmony, 
                              prefix = prefix, 
                              vars.to.regress = vars.to.regress,
                              bubblemap_file = bubblemap_file,
                              essential_genes = essential_genes)
  obj=reslist1$obj
  files=reslist1$files
  rm(reslist1)

  write.csv(files, paste0(prefix, ".iter_png.csv"))

  setwd(cur_folder)

  celltypes<-unique(obj$layer4)
  celltypes<-celltypes[order(celltypes)]
  ctdf<-data.frame("celltype"=celltypes, "resolution"=0.01)
  write.table(ctdf, paste0(prefix, ".scDynamic.celltype_res.txt"), row.names=F, sep="\t", quote=F)

  obj<-factorize_layer(obj, "layer4")
  Idents(obj)<-"layer4"

  saveRDS(obj@meta.data, paste0(prefix, ".scDynamic.meta.rds"))
  write.csv(obj@meta.data, paste0(prefix, ".scDynamic.meta.csv"))

  save_umap(paste0(prefix, ".scDynamic.umap"), obj)

  if(length(unique(obj$layer4)) > 1){
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

    width<-max(5000, min(10000, length(unique(Idents(obj))) * 400 + 1000))
    height<-max(3000, min(10000, length(all_top10) * 60 + 1000))
    png(paste0(prefix, ".layer4.heatmap.png"), width=width, height=height, res=300)
    print(g)
    dev.off()
  }
  
  output_celltype_figures(obj, "layer4", prefix, bubblemap_file, cell_activity_database, combined_ct_source, group.by="orig.ident", name="sample")

  if(!by_individual_sample){
    #output individual sample dot plot, with global scaled average gene expression.
    obj$sample_layer4<-paste0(obj$orig.ident, ":", obj$layer4)
    g<-get_bubble_plot( obj = obj, 
                    cur_res = NA, 
                    cur_celltype = "sample_layer4", 
                    bubblemap_file = bubblemap_file, 
                    assay = "RNA", 
                    orderby_cluster = F)
    gdata<-g$data
    gdata$id<-gsub(".+: ","",gdata$id)
    gdata$sample<-gsub(":.+","",gdata$id)
    gdata$id<-gsub(".+:","",gdata$id)

    for(sample in unique(gdata$sample)){
      sdata<-gdata[gdata$sample == sample,]
      g$data=sdata
      dot_file = paste0(prefix, ".layer4.", sample, ".dot.png")
      png(dot_file, width=4000, height=get_dot_height(obj, "layer4"), res=300)
      print(g)
      dev.off()
    }

    has_batch<-FALSE
    if("batch" %in% colnames(obj@meta.data)){
      if("sample" %in% colnames(obj@meta.data)){
        has_batch=any(obj$batch != obj$sample)
      }else{
        has_batch=any(obj$batch != obj$orig.ident)
      }
    }
    if(has_batch){
      output_celltype_figures(obj, "layer4", prefix, bubblemap_file, cell_activity_database, combined_ct_source, group.by="batch", name="batch")
    }
  }

  obj$seurat_layer4=paste0(obj$layer4_clusters, ": ", obj$layer4_raw)

  cur_celltype="layer4"
  for(pct in unique(unlist(obj[[cur_celltype]]))){
    cells=colnames(obj)[obj[[cur_celltype]] == pct]
    subobj=subset(obj, cells=cells)
    subobj$seurat_layer4=paste0(subobj$layer4_clusters, ": ", subobj$layer4_raw)
    g<-get_dim_plot(subobj, group.by="layer4_clusters", label.by="seurat_layer4", reduction="umap", legend.title="")

    png(paste0(prefix, ".", cur_celltype, ".", celltype_to_filename(pct), ".umap.png"), width=2400, height=2000, res=300)
    print(g)
    dev.off()
  }

  # obj$ct_cluster<-paste0(obj$layer4, ":", obj$layer4_clusters)
  # dot_width=4400
  # g<-get_bubble_plot(obj, cur_res="layer4_clusters", "layer4", bubblemap_file, assay="RNA", orderby_cluster = FALSE)
  # dot_file = paste0(prefix, ".layer4.dot.sub.png")
  # png(dot_file, width=dot_width, height=get_dot_height(obj, "ct_cluster"), res=300)
  # print(g)
  # dev.off()

  tb=data.frame("cell"=colnames(obj), "cell_type"=obj$layer4, category=prefix)

  output_file=paste0(prefix,".dynamic.html")
  rmdfile = "seurat_scDynamic_one_layer_one_resolution.rmd"
  rmarkdown::render(rmdfile, output_file=output_file)
  return(list(html=output_file, ct_count=tb))
}

if(by_individual_sample){
  if("harmony" %in% names(obj@reductions)){
    obj@reductions["harmony"]<-NULL
  }
  if("umap" %in% names(obj@reductions)){
    obj@reductions["umap"]<-NULL
  }
  if("pca" %in% names(obj@reductions)){
    obj@reductions["pca"]<-NULL
  }
  root_folder = getwd()
  samples <- unique(obj$orig.ident)
  result_list = c()
  all_ct_counts = NULL

  sample=samples[2]
  for (sample in samples){
    cat("processing ", sample, "...\n")
    cur_folder = paste0(root_folder, "/", sample)
    if(!dir.exists(cur_folder)){
      dir.create(cur_folder)
    }
    tmp_folder = paste0(cur_folder, "/details")
    if(!dir.exists(tmp_folder)){
      dir.create(tmp_folder)
    }
    
    f1<-read.table(paste0(root_folder, "/fileList1.txt"), sep="\t")
    f1$V1[f1$V2=="task_name"] = sample
    write.table(f1, paste0(cur_folder, "/fileList1.txt"), sep="\t", row.names=F)

    file.copy(paste0(root_folder, "/scRNA_func.r"), paste0(cur_folder, "/scRNA_func.r"), overwrite = TRUE)
    file.copy(paste0(root_folder, "/seurat_scDynamic_one_layer_one_resolution.rmd"), paste0(cur_folder, "/seurat_scDynamic_one_layer_one_resolution.rmd"), overwrite = TRUE)
    file.copy(paste0(root_folder, "/reportFunctions.R"), paste0(cur_folder, "/reportFunctions.R"), overwrite = TRUE)

    subobj<-subset(obj, cells=colnames(obj)[obj@meta.data$orig.ident==sample])

    #remove genes without count
    sub_counts<-subobj@assays$RNA@counts
    rc<-rowSums(sub_counts)
    rc<-rc[rc>0]
    rm(sub_counts)
    subobj <- subset(subobj, features = names(rc))

    if(0){#for test
      subobj=subset(subobj, cells=colnames(subobj)[1:2000])
    }

    #for each sample, do its own PCA, FindClusters and UMAP first
    subobj = sub_cluster( subobj = subobj, 
                          assay = assay, 
                          by_sctransform = by_sctransform, 
                          by_harmony = FALSE, 
                          redo_harmony = FALSE, 
                          curreduction = "pca", 
                          k_n_neighbors = min(npcs, 20),
                          u_n_neighbors = min(npcs, 30),
                          random.seed = random.seed,
                          resolutions = resolution,
                          cur_npcs = npcs, 
                          cur_pca_dims = c(1:npcs),
                          vars.to.regress = vars.to.regress,
                          essential_genes = essential_genes)    

    res_list = do_analysis( tmp_folder = tmp_folder,
                            cur_folder = cur_folder,
                            obj = subobj, 
                            layer2map = layer2map, 
                            npcs = npcs, 
                            resolution = resolution, 
                            random.seed = random.seed, 
                            by_sctransform = by_sctransform, 
                            by_harmony = 0, 
                            prefix = sample, 
                            vars.to.regress = vars.to.regress, 
                            bubblemap_file = bubblemap_file, 
                            essential_genes = essential_genes,
                            by_individual_sample = 1)

    result_list<-c(result_list, res_list$html)
    all_ct_counts<-rbind(all_ct_counts, res_list$ct_count)
  }
  setwd(root_folder)
  writeLines(result_list, paste0(prefix, ".samples.list"))
  write.table(all_ct_counts, paste0(prefix, ".ct_count.tsv"), sep="\t", row.names=T)

  output_file=paste0(prefix,".dynamic_individual.html")
  rmdfile = "seurat_scDynamic_one_layer_one_resolution_summary.rmd"
  rmarkdown::render(rmdfile, output_file=output_file)
}else{
  cur_folder = getwd()
  tmp_folder = paste0(cur_folder, "/details")
  if(!dir.exists(tmp_folder)){
    dir.create(tmp_folder)
  }
  res_list <- do_analysis(tmp_folder = tmp_folder,
                          cur_folder = cur_folder,
                          obj = obj, 
                          layer2map = layer2map, 
                          npcs = npcs, 
                          resolution = resolution, 
                          random.seed = random.seed, 
                          by_sctransform = by_sctransform, 
                          by_harmony = by_harmony, 
                          prefix = prefix, 
                          vars.to.regress = vars.to.regress, 
                          bubblemap_file = bubblemap_file, 
                          essential_genes = essential_genes,
                          by_individual_sample = 0);
}
