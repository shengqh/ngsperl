source("scRNA_func.r")

library(dplyr)
library(Seurat)
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(rmdformats)
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

bubblemap_file=myoptions$bubblemap_file

if(file.exists(parFile2)){
  npcs<-read.table(parFile2, row.names=1)$V2[1]
}
pca_dims<-1:npcs

tiers<-read.table(myoptions$HLA_panglao5_file, sep="\t", header=T)

remove_subtype_of=remove_subtype
cell_activity_database<-read_cell_markers_file(markerfile, species, remove_subtype_of, HLA_panglao5_file, curated_markers_file=myoptions$curated_markers_file)

prefix<-outFile

if(!exists("obj")){
  finalList<-readRDS(parFile1)
  obj<-finalList$obj
}

resolutions=seq(from = 0.1, to = 0.9, by = 0.1)
#resolutions=seq(from = 0.1, to = 0.3, by = 0.1)

layer1map<-split(tiers$Layer1, tiers$Celltype.name)
layer2map<-split(tiers$Layer2, tiers$Celltype.name)
layer3map<-split(tiers$Layer3, tiers$Celltype.name)
layer4map<-split(tiers$Layer4, tiers$Celltype.name)

obj[["layer0"]]<-"Unassigned"

previous_layer<-"layer0"
cur_layer="layer1"
layermap=layer1map

layer_cluster_celltype<-function(obj, previous_layer, cur_layer, layermap, npcs, resolutions, random.seed, by_sctransform, by_harmony){
  meta<-obj@meta.data
  previous_celltypes<-unique(meta[[previous_layer]])
  
  assay=ifelse(by_sctransform, "SCT", "RNA")
  clusters=paste0(assay, "_snn_res.", resolutions)
  
  allcts<-NULL
  pct<-previous_celltypes[length(previous_celltypes)]
  for(pct in previous_celltypes){
    cat("celltype", pct, "...\n")
    cells<-rownames(meta)[meta[[previous_layer]] == pct]
    subobj<-subset(obj, cells=cells)
    
    if(length(cells) <= 50) {
      #assign previous cell type
      cts<-subobj[[previous_layer]]
      oldct<-subobj@meta.data[,previous_layer]
      for(cluster in clusters){
        cts[,cluster]<-oldct
      }
      cts$layer<-oldct
    }else{
      if(by_harmony){
        cat("harmony\n")
        DefaultAssay(subobj)<-"RNA"
        subobj[['harmony']]<-NULL
        subobj<-do_harmony(subobj, npcs, batch_file="")
        curreduction="harmony"
      }else if (by_sctransform) {
        cat("sctransform\n")
        DefaultAssay(subobj)<-"RNA"
        subobj[['SCT']]<-NULL
        subobj<-SCTransform(subobj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
        subobj<-RunPCA(subobj)
        curreduction="pca"
      }else{
        cat("normalization\n")
        DefaultAssay(subobj)<-"RNA"
        subobj<-NormalizeData(subobj)
        subobj<-FindVariableFeatures(subobj)
        subobj<-ScaleData(subobj, vars.to.regress = c("percent.mt"))
        subobj<-RunPCA(subobj)
        curreduction="pca"
      }
      cat("FindClusters\n")
      subobj<-FindNeighbors(object=subobj, reduction=curreduction, dims=pca_dims, verbose=FALSE)
      subobj<-FindClusters(object=subobj, random.seed=random.seed, resolution=resolutions, verbose=FALSE)
      
      cat("RunUMAP\n")
      subobj<-RunUMAP(object = subobj, reduction=curreduction, dims=pca_dims, verbose = FALSE)
      
      cat("Cell type annotation\n")
      cts<-subobj[[previous_layer]]
      cluster = clusters[1]
      for(cluster in clusters){
        cat("  ", cluster, "\n")
        data.norm=get_seurat_average_expression(subobj, cluster)
        
        predict_celltype<-ORA_celltype(data.norm,cell_activity_database$cellType,cell_activity_database$weight)
        
        new.cluster.ids<-names(predict_celltype$max_cta)
        #names(new.cluster.ids) <- colnames(data.norm)
        
        layer_ids<-unlist(layermap[new.cluster.ids])
        names(layer_ids) <- colnames(data.norm)
        
        oldcluster<-subobj[[cluster]][[1]]
        newct<-layer_ids[oldcluster]
        newdf<-data.frame("ct"=newct, "cl"=oldcluster)
        newdf$ct_cl<-paste0(newdf$ct, ": ", newdf$cl)
        colnames(newdf)<-paste0(cluster, "_", colnames(newdf))
        cts<-cbind(cts, newdf)
      }
      
      cls<-cts[,grepl(paste0(assay, "_snn_res.+_ct$"), colnames(cts))]
      
      votes<-apply(cls, 1, function(x){
        tb<-table(unlist(x))
        tbwhich<-which(tb==max(tb))
        return(names(tb)[tbwhich[1]])
      })
      cts$layer<-unlist(votes)
      subobj$layer=cts$layer

      # votes2<-apply(cls, 1, function(x){
      #   df<-data.frame("ct"=unlist(x), "weight"=resolutions)
      #   adf<-aggregate(df$weight, by=list(celltype=df$ct), FUN=sum)
      #   adf<-adf[order(adf$x, decreasing = T),]
      #   return(adf$celltype[1])
      # })
      # cts$layer<-unlist(votes2)
      # subobj$layer=cts$layer

      g<-DimPlot(subobj, group.by = "layer", label=T) + ggtitle(pct)
      if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
        g2<-get_bubble_plot(subobj, NA, "layer", bubblemap_file)
        g<-g+g2+plot_layout(ncol = 2, widths = c(3, 5))
        width=6400
      }else{
        width=2400
      }
      png(paste0(prefix, ".", previous_layer, "_", gsub(" ", "_", pct), ".png"), width=width, height=2000, res=300)
      print(g)
      dev.off()
      
      g<-NULL
      for (cl in colnames(cls)){
        ct<-paste0(cl, "_cl")
        subobj[[ct]]=unlist(cts[,ct])
        gg<-DimPlot(subobj, group.by = ct, label=T) + ggtitle(cl)
        if(is.null(g)){
          g=gg
        }else{
          g<-g+gg
        }
      }
      png(paste0(prefix, ".", previous_layer, "_", gsub(" ", "_", pct), ".detail.png"), width=7000, height=4800, res=300)
      g<-g+plot_layout(ncol=3)
      print(g)
      dev.off()
    }
    
    allcts<-rbind(allcts, cts)
    
    if(previous_layer == "layer0"){
      obj[['umap']] = subobj[['umap']]
    }
    rm(subobj)
  }
  obj[[cur_layer]]=unlist(allcts[colnames(obj), "layer"])
  write.csv(allcts, file=paste0(previous_layer, "_", cur_layer, ".csv"))
  
  g<-DimPlot(obj, group.by = cur_layer, label=T)
  if(!is.null(bubblemap_file) && file.exists(bubblemap_file)){
    g2<-get_bubble_plot(obj, NA, cur_layer, bubblemap_file)
    g<-g+g2+plot_layout(ncol = 2, widths = c(3, 5))
    width=6400
  }else{
    width=2400
  }
  png(paste0(prefix, ".", previous_layer, "_to_", cur_layer, ".png"), width=width, height=2000, res=300)
  print(g)
  dev.off()
  
  return(obj)
}

obj<-layer_cluster_celltype(obj, "layer0", "layer1", layer1map, npcs, resolutions, random.seed, by_sctransform, FALSE)
obj<-layer_cluster_celltype(obj, "layer1", "layer2", layer2map, npcs, resolutions, random.seed, by_sctransform, FALSE)
obj<-layer_cluster_celltype(obj, "layer2", "layer3", layer3map, npcs, resolutions, random.seed, by_sctransform, FALSE)
obj<-layer_cluster_celltype(obj, "layer3", "layer4", layer4map, npcs, resolutions, random.seed, by_sctransform, FALSE)

write.csv(obj@meta.data, paste0(prefix, ".meta.csv"))

saveRDS(obj, paste0(prefix, ".multires.rds"))
