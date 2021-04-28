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
require(data.table)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

optionstr<-paste0(options_table$V1, collapse = "_")

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
Remove_Mt_rRNA= ifelse(myoptions$Remove_Mt_rRNA == "FALSE", FALSE, TRUE)
resolution=as.numeric(myoptions$resolution)
pool_sample<-ifelse(myoptions$pool_sample == "0", FALSE, TRUE)
batch_for_integration<-ifelse(myoptions$batch_for_integration == "0", FALSE, TRUE)
by_integration<-ifelse(myoptions$by_integration == "0", FALSE, TRUE)
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)

hto_str = ""
hto_data = list()
hto_sample_file<-myoptions$hto_sample_file
has_hto = (hto_sample_file != "")
if (has_hto) {
  hto_samples = read.table(hto_sample_file, sep="\t", header=T, stringsAsFactors=FALSE)
  hto_md5 = list()
  hto_md5[["hto_samples"]] = digest(hto_sample_file, file=TRUE)
  if (!exists('parSampleFile5')) {
    stop('hto_cell_file is not defined.')
  }
  if (!file.exists(parSampleFile5)) {
    stop(paste0('hto_cell_file ', parSampleFile5, ' is not defined or not exists.'))
  }
  hto_cell_files = read.table(parSampleFile5, sep="\t", stringsAsFactors=FALSE, row.names=2)
  #sample = rownames(hto_cell_files)[1]
  for (sample in rownames(hto_cell_files)){
    cell_file = hto_cell_files[sample, "V1"]
    cell_data = read.csv(cell_file, stringsAsFactors=FALSE, header=TRUE, check.names=F)
    hto_md5[[sample]] = digest(cell_file, file=TRUE)
    cell_data$Sample = ""
    cur_samples = hto_samples[hto_samples$File == sample,]
    for (cidx in c(1:nrow(cur_samples))){
      tagname = cur_samples$Tagname[cidx]
      samplename = cur_samples$Sample[cidx]
      cell_data$Sample[cell_data$HTO == tagname] = samplename
    }
    sample_cell_data = cell_data[cell_data$Sample != "",]
    hto_data[[sample]] = sample_cell_data
  }
  hto_str = paste0(hto_md5, collapse="_")
}
prefix<-outFile

nFeature_cutoff_min=as.numeric(myoptions$nFeature_cutoff_min)
nFeature_cutoff_max=as.numeric(myoptions$nFeature_cutoff_max)
nCount_cutoff=as.numeric(myoptions$nCount_cutoff)
mt_cutoff=as.numeric(myoptions$mt_cutoff)
species=myoptions$species
#nCount_sd_cutoff=as.numeric(options$nCount_sd_cutoff)
nCount_sd_cutoff=0

pca_dims<-1:as.numeric(myoptions$pca_dims)

filelist1<-read.table(parSampleFile1, header=F)
filestr<-paste0(filelist1$V1, collapse = "_")

prefixMd5=digest(paste0(optionstr, filestr, hto_str), "md5", serialize = FALSE)

finalList<-list()
finalListFile<-paste0(prefix, ".final.rds")
finalListFileOptionsMd5<-paste0(prefix, ".final.options.md5")

isCalc=TRUE
if(file.exists(finalListFile) && file.exists(finalListFileOptionsMd5)){
  oldMd5=str_trim(readChar(finalListFileOptionsMd5, file.info(finalListFileOptionsMd5)$size))
  isCalc=(oldMd5 != prefixMd5)
}

if(isCalc){
  #read raw count dat
  filelist1<-read.table(parSampleFile1, header=F, stringsAsFactors = F)
  rawobjs = list()
  fidx=1
  for (fidx in c(1:nrow(filelist1))) {
    fileName  = filelist1[fidx, 1]
    fileTitle = filelist1[fidx, 2]
    if(dir.exists(fileName)){
      counts = Read10X(fileName)
    } else {
      counts = Read10X_h5(fileName)
    }
    
    adt.counts<-NULL
    if (is.list(counts)){
      adt.counts<-counts$`Antibody Capture`
      counts<-counts$`Gene Expression` 
    }

    if (species=="Mm") {
      rownames(counts)<-toMouseGeneSymbol(rownames(counts))
    }
    if (species=="Hs") {
      rownames(counts)<-toupper(rownames(counts))
    }
    sobj = CreateSeuratObject(counts = counts, project = fileTitle)
    sobj[["percent.mt"]] <- PercentageFeatureSet(object = sobj, pattern = Mtpattern)
    if (!is.null(adt.counts)){
      mat<-as.matrix(adt.counts)
      rowsum<-apply(mat>0, 1, sum)
      mat<-mat[rowsum > (ncol(mat) / 2),]
      sobj[["ADT"]] <- CreateAssayObject(counts = mat)
    }

    if(has_hto && fileTitle %in% names(hto_data)) {
      cell_data = hto_data[[fileTitle]]
      validobj = subset(sobj, cells=cell_data$X)
      #tagname = unique(cell_data$HTO)[1]
      for (tagname in unique(cell_data$HTO)){
        tagcells = cell_data[cell_data$HTO == tagname,]
        sample = tagcells$Sample[1]
        sample_obj = subset(validobj, cells=tagcells$X)
        sample_obj$orig.ident = sample
        sample_obj<-RenameCells(object=sample_obj, new.names=paste0(sample, "_", colnames(sample_obj)))
        rawobjs[[sample]] = sample_obj
      }
    }else{
      sobj<-RenameCells(object=sobj, new.names=paste0(fileTitle, "_", colnames(sobj)))
      rawobjs[[fileTitle]] = sobj
    }
  }
  
  if(pool_sample){
    if (!exists('parSampleFile4')) {
      stop('pool_sample_groups is not defined.')
    }
    if (!file.exists(parSampleFile4)) {
      stop(paste0('pool_sample_groups ', parSampleFile4, ' is not defined or not exists.'))
    }
    #cat("pooling samples ... \n")
    pools = read.table(parSampleFile4, header=F, stringsAsFactors = F)
    poolNames = unique(pools$V2)
    pooledObjs = lapply(poolNames, function(pn){
      curPools<-pools[pools$V2==pn,]
      curListIndecies<-which(unlist(lapply(rawobjs, function(x) x@project.name %in% curPools$V1)))
      curObjs<-rawobjs[curListIndecies]
      if(length(curObjs) == 1){
        curobj = curObjs[[1]]
      }else{
        curobj <- merge(curObjs[[1]], y = unlist(curObjs[2:length(curObjs)]), project = "integrated")
      }
      curobj@project.name=pn
      curobj$orig.ident=rep(pn, length(curobj$orig.ident))
      Idents(curobj)<-"orig.ident"
      return(curobj)
    });
    notPoolIndecies<-which(unlist(lapply(rawobjs, function(x) !(x@project.name %in% pools$V1))))
    if(length(notPoolIndecies) > 0){
      notPoolObjs<-rawobjs[notPoolIndecies]
      rawobjs<-c(pooledObjs, notPoolObjs)
    }else{
      rawobjs<-pooledObjs
    }
    #cat("pooling samples done. \n")
  }
  
  if(length(rawobjs) == 1){
    rawobj <- rawobjs[[1]]
  }else{
    rawobj <- merge(rawobjs[[1]], y = unlist(rawobjs[2:length(rawobjs)]), project = "integrated")
  }

  finalList$rawobj<-rawobj

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

  if(by_integration){
    if(batch_for_integration){
      if(!file.exists(parSampleFile3)){
        stop("batch_for_integration_groups not defined.")
      }
      pools = read.table(parSampleFile3, header=F, stringsAsFactors = F)
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

    if(by_sctransform){
      #perform sctransform
      objs<-lapply(objs, function(x){
        x <- SCTransform(x, verbose = FALSE)
        return(x)
      })  
      
      if(length(objs) == 1){
        obj <- objs[[1]]
      }else{
        all_genes <- Reduce(intersect, lapply(objs, rownames)) 
        obj_features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
        objs <- PrepSCTIntegration(object.list = objs, anchor.features = obj_features, verbose = FALSE)
        obj_anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT", anchor.features = obj_features, verbose = FALSE)
        obj <- IntegrateData(anchorset = obj_anchors, normalization.method = "SCT", verbose = FALSE, features.to.integrate=all_genes)      
      }
    }else{
      #perform standard workflow
      objs<-lapply(objs, function(x){
        x <- NormalizeData(x, verbose = FALSE)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)  
        return(x)
      })  
      
      if(length(objs) == 1){
        obj <- objs[[1]]
      }else{
        all_genes <- Reduce(intersect, lapply(objs, rownames)) 
        obj.anchors <- FindIntegrationAnchors(object.list = objs, dims = 1:20)
        obj <- IntegrateData(anchorset = obj.anchors, dims = 1:20, features.to.integrate=all_genes)    
        obj <- ScaleData(obj, verbose = FALSE)
      }
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
  }
  rm(objs)
  
  obj_markers <- run_cluster(obj, Remove_Mt_rRNA, rRNApattern, Mtpattern, pca_dims, by_sctransform)

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
}
