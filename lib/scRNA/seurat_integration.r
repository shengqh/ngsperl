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
library(patchwork)
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

prefix<-outFile

has_batch_file<-file.exists(parSampleFile2)

pca_dims<-1:as.numeric(myoptions$pca_dims)

finalList<-list()
finalListFile<-paste0(prefix, ".preprocess.rds")

rawobj<-readRDS(parFile1)

finalList<-preprocessing_rawobj(rawobj, myoptions, prefix)
rawobj<-finalList$rawobj
finalList<-finalList[names(finalList) != "rawobj"]

poolmap = get_batch_samples(parSampleFile2, unique(rawobj$sample))

rawobj$batch=unlist(poolmap[rawobj$sample])
objs<-SplitObject(object = rawobj, split.by = "batch")
rm(rawobj)

#filter cells
finalList$filter<-list(nFeature_cutoff_min=nFeature_cutoff_min,
                      nFeature_cutoff_max=nFeature_cutoff_max,
                      mt_cutoff=mt_cutoff,
                      nCount_cutoff=nCount_cutoff,
                      nCount_sd_cutoff=nCount_sd_cutoff)

if(by_sctransform){
  cat("performing SCTransform ...")
  #perform sctransform
  objs<-lapply(objs, function(x){
    x <- SCTransform(x, verbose = FALSE)
    return(x)
  })  
  obj_features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
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

if(length(objs) == 1){
  obj <- objs[[1]]
}else{
  if(by_sctransform){
    all_genes <- Reduce(intersect, lapply(objs, rownames)) 
    objs <- PrepSCTIntegration(object.list = objs, anchor.features = obj_features, verbose = FALSE)
    obj_anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT", anchor.features = obj_features, verbose = FALSE)
    obj <- IntegrateData(anchorset = obj_anchors, normalization.method = "SCT", verbose = FALSE, features.to.integrate=all_genes)      
    obj <- RunPCA(object = obj, verbose=FALSE)
  }else{
    all_genes <- Reduce(intersect, lapply(objs, rownames)) 
    obj.anchors <- FindIntegrationAnchors(object.list = objs, dims = 1:20)
    obj <- IntegrateData(anchorset = obj.anchors, dims = 1:20, features.to.integrate=all_genes)  
    DefaultAssay(obj) <- "integrated"  
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(object = obj, verbose=FALSE)
  }
}

obj <- run_cluster_only(obj, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform, reduction="pca")

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

output_integraion_dimplot(obj, outFile, has_batch_file)
