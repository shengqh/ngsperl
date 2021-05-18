library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(cowplot)
library(scales)
library(stringr)
library(harmony)
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

optionstr<-paste0(options_table$V1, collapse = "_")

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_Mt_rRNA= ifelse(myoptions$Remove_Mt_rRNA == "FALSE", FALSE, TRUE)
resolution=as.numeric(myoptions$resolution)
batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_integration<-ifelse(myoptions$by_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
integration_by_harmony=ifelse(myoptions$integration_by_harmony == "0", FALSE, TRUE)

prefix<-outFile

nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
mt_cutoff=as.numeric(myoptions$mt_cutoff)
species=myoptions$species
#nCount_sd_cutoff=as.numeric(options$nCount_sd_cutoff)
nCount_sd_cutoff=0

pca_dims<-1:as.numeric(myoptions$pca_dims)

finalList<-list()
finalListFile<-paste0(prefix, ".preprocess.rds")

rawobjs<-readRDS(parFile1)

#filter cells
finalList$filter<-list(nFeature_cutoff_min=nFeature_cutoff_min,
                      nFeature_cutoff_max=nFeature_cutoff_max,
                      mt_cutoff=mt_cutoff,
                      nCount_cutoff=nCount_cutoff,
                      nCount_sd_cutoff=nCount_sd_cutoff)

objs<-lapply(rawobjs, function(x){
  sobj<-subset(x, subset = nFeature_RNA > nFeature_cutoff_min & nFeature_RNA<nFeature_cutoff_max & nCount_RNA > nCount_cutoff & percent.mt < mt_cutoff)

  if(nCount_sd_cutoff > 0){
    nCount_mean = mean(sobj[["nCount_RNA"]])
    nCount_sd = sd(sobj[["nCount_RNA"]])
    nCount_sd_min = nCount_mean - nCount_sd_cutoff * nCount_sd
    nCount_sd_max = nCount_mean + nCount_sd_cutoff * nCount_sd
    sobj<-subset(sobj, subset = nCount_RNA > nCount_sd_min & nCount_RNA < nCount_sd_max)
    finalList$filter$nCount_sd_min = nCount_sd_min
    finalList$filter$nCount_sd_max = nCount_sd_max
  }
  return(sobj)
})  
rm(rawobjs)

if(by_sctransform){
  cat("performing SCTransform ...")
  #perform sctransform
  objs<-lapply(objs, function(x){
    x <- SCTransform(x, verbose = FALSE)
    return(x)
  })  
  assay="SCT"
}else{
  cat("performing NormalizeData/FindVariableFeatures ...")
  #perform standard workflow
  objs<-lapply(objs, function(x){
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)  
    return(x)
  })  
  assay="RNA"
}

integration
if(length(objs) == 1){
  obj=obj[[1]]
}else{
  if(by_integration){
    if(integration_by_harmony){
      if(file.exists(parSampleFile2)){
        pools = read.table(parSampleFile2, header=F, stringsAsFactors = F)
        poolmap=split(pools$V2, pools$V1)
      }else{
        poolmap=split(names(objs),names(objs))
      }
      
      if(by_sctransform){
        objs.features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)  
        obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
        VariableFeatures(obj) <- objs.features        
      }else{
        obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
        obj <- NormalizeData(obj, verbose = FALSE)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
        all.genes <- rownames(obj)  
        obj <- ScaleData(obj, features = all.genes, verbose = FALSE)
      }
      obj <- RunPCA(object = obj, assay=assay, verbose=FALSE)
      obj$batch <- unlist(poolmap[obj$orig.ident])
      obj <- RunHarmony(object = obj,
                        assay.use = "SCT",
                        reduction = "pca",
                        dims.use = pca_dims,
                        group.by.vars = "batch")
      reduction="harmony"
    }else{
      if(file.exists(parSampleFile2)){
        pools = read.table(parSampleFile2, header=F, stringsAsFactors = F)
        pools$V2=as.character(pools$V2)
        poolNames = unique(pools$V2)
        pooledObjs = lapply(poolNames, function(pn){
          curPools<-pools[pools$V2==pn,]
          curListIndecies<-which(unlist(lapply(objs, function(x) x@project.name %in% curPools$V1)))
          curObjs<-objs[curListIndecies]
          if(length(curObjs) ==1){
            curobj<-curObjs[[1]]
          }else{
            curobj <- merge(curObjs[[1]], y = unlist(curObjs[2:length(curObjs)]), project = "integrated")
          }
          curobj@project.name=pn
          return(curobj)
        });
        notPoolIndecies<-which(unlist(lapply(objs, function(x) !(x@project.name %in% pools$V1))))
        if(length(notPoolIndecies) > 0){
          notPoolObjs<-objs[notPoolIndecies]
          objs<-c(pooledObjs, notPoolObjs)
        }else{
          objs<-pooledObjs
        }
      }
  
      if(length(objs) == 1){
        obj <- objs[[1]]
      }else{
        if(by_sctransform){
          all_genes <- Reduce(intersect, lapply(objs, rownames)) 
          obj_features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
          objs <- PrepSCTIntegration(object.list = objs, anchor.features = obj_features, verbose = FALSE)
          obj_anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT", anchor.features = obj_features, verbose = FALSE)
          obj <- IntegrateData(anchorset = obj_anchors, normalization.method = "SCT", verbose = FALSE, features.to.integrate=all_genes)      
          obj <- RunPCA(object = obj, verbose=FALSE)
        }else{
          all_genes <- Reduce(intersect, lapply(objs, rownames)) 
          obj.anchors <- FindIntegrationAnchors(object.list = objs, dims = 1:20)
          obj <- IntegrateData(anchorset = obj.anchors, dims = 1:20, features.to.integrate=all_genes)    
          obj <- ScaleData(obj, verbose = FALSE)
          var.genes <- VariableFeatures(obj)
          obj <- RunPCA(object = obj, features = var.genes, verbose=FALSE)
        }
      }
      reduction="pca"
    }
  }else{
    if(length(objs) == 1){
      obj <- objs[[1]]
    }else{
      obj <- merge(objs[[1]], y = unlist(objs[2:length(objs)]), project = "integrated")
    }
    if(by_sctransform){
      #perform sctransform
      obj <- SCTransform(obj, verbose = FALSE)
    }else{
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      all.genes <- rownames(obj)  
      obj <- ScaleData(obj, features = all.genes, verbose = FALSE)
    }
    reduction="pca"
  }
}

rm(objs)

obj_markers <- run_cluster(obj, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, reduction=reduction)

clusterMarkers<-obj_markers$markers %>% group_by(cluster)
write.csv(clusterMarkers, file=paste0(prefix, ".allmarkers.csv"), row.names=F, quote = F)

finalList$markers<-obj_markers$markers

obj<-obj_markers$object
seurat_clusters<-unlist(obj[["seurat_clusters"]])
seurat_colors<-hue_pal()(length(levels(seurat_clusters)))
names(seurat_colors)<-levels(seurat_clusters)

finalList$seurat_colors<-seurat_colors

counts=GetAssayData(obj,assay="RNA",slot="counts")
clusters=obj@active.ident
sumcounts=get_cluster_count(counts, clusters)
logsumcounts<-log2(sumcounts+1)

data.quantileAll<-apply(logsumcounts, 2, function(x){quantile(x, 0.75)})

norm_method=""
if(any(data.quantileAll == 0)){
  norm_method = ".normByTotal"
  data.all <- apply(logsumcounts, 2, sum)
  data.all<-data.all / median(data.all)
  data.norm <- t(t(logsumcounts) / data.all)
}else{
  norm_method = ".normByUpQuantile"
  data.quantileAll<-data.quantileAll / median(data.quantileAll)
  data.norm <- t(t(logsumcounts) / data.quantileAll)
}

write.csv(sumcounts, file=paste0(prefix, ".cluster.count.csv"))
write.csv(data.norm, file=paste0(prefix, ".cluster", norm_method, ".csv"))

clusters<-data.frame("cell" = c(1:length(obj$seurat_clusters)), "seurat_clusters"=as.numeric(as.character(obj$seurat_clusters)), stringsAsFactors = F)
rownames(clusters)<-names(obj$seurat_clusters)
write.csv(clusters, file=paste0(prefix, ".cluster.csv"))

finalList$obj<-obj
saveRDS(finalList, file=finalListFile)

fileConn<-file(finalListFileOptionsMd5)
writeLines(c(prefixMd5), fileConn)
close(fileConn)
