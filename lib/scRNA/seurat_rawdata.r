rm(list=ls()) 
outFile='AG_integrated'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20220803_integrated_project/seurat_rawdata/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(digest)
library(patchwork)
library(sparseMatrixStats)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
hemoglobinPattern=myoptions$hemoglobinPattern

species=myoptions$species
pool_sample<-ifelse(myoptions$pool_sample == "0", FALSE, TRUE)

hto_str = ""
hto_data = list()
hto_sample_file<-myoptions$hto_sample_file
has_hto = (hto_sample_file != "")

if (has_hto) {
  hto_samples = read.table(hto_sample_file, sep="\t", header=T, stringsAsFactors=FALSE)
  hto_md5 = list()
  hto_md5[["hto_samples"]] = digest(hto_sample_file, file=TRUE)
  if (!exists('parSampleFile4')) {
    stop('hto_cell_file is not defined.')
  }
  if (!file.exists(parSampleFile4)) {
    stop(paste0('hto_cell_file ', parSampleFile4, ' is not defined or not exists.'))
  }
  hto_cell_files = read.table(parSampleFile4, sep="\t", stringsAsFactors=FALSE, row.names=2)
  sample = rownames(hto_cell_files)[1]
  for (sample in rownames(hto_cell_files)){
    cell_file = hto_cell_files[sample, "V1"]
    if(file_ext(cell_file) == "csv"){
      cell_data = read.csv(cell_file, stringsAsFactors=FALSE, header=TRUE, check.names=F)
      colnames(cell_data)[1]<-"cell"
    }else{
      cell_data = readRDS(cell_file)
      cell_data$cell<-rownames(cell_data)
      cell_data$HTO<-cell_data$final
    }
    cell_data$Sample = ""
    cur_samples = hto_samples[hto_samples$File == sample,]
    if(!all(cur_samples$Tagname %in% cell_data$HTO)){
      stop(paste0("failed, not all tag defined ", paste0(cur_samples$Tagname, collapse = "/") , " for sample ", sample,  " were found in HTO tags ", paste0(unique(cell_data$HTO), collapse = "/") ))
    }
    
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

remove_cells<-list()
if(file.exists(parFile1)){
  load(parFile1)
  
  objs_config=read.table(parFile2, sep="\t", header=T, stringsAsFactors = F)
  
  if(is.null(names(object.list))){
    names(object.list)<-objs_config$sample
  }
  
  rowi=1
  for(rowi in c(1:nrow(objs_config))){
    sample_id = objs_config$sample[rowi]

    #cat(sample_id, "\n")
    
    if(sample_id %in% names(object.list)){
      crs<-objs_config$cluster_remove[rowi]
      if(!is.na(crs)){
        sample_obj = object.list[[sample_id]]
        cd<-sample_obj$meta
  
        crslist<-unlist(strsplit(crs,','))
        rcs<-rownames(cd)[as.character(cd$seurat_clusters) %in% crslist] 
        
        remove_cells[[sample_id]]=rcs
      }
    }
  }
}

read_gzip_count_file<-function(files, sample, species){
  all_counts<-NA
  index = 1
  for(file in files){
    cat('reading ', file, "\n")
    counts<-fread(file)
    if (species=="Mm") {
      counts$GENE<-toMouseGeneSymbol(counts$GENE)
    }
    if (species=="Hs") {
      counts$GENE<-toupper(counts$GENE)
    }
    
    counts<-counts[!duplicated(counts$GENE),]
    colnames(counts)[2:ncol(counts)]<-paste0(sample, index, "_", colnames(counts)[2:ncol(counts)])
    if(is.na(all_counts)){
      all_counts<-counts
    }else{
      all_counts<-merge(all_counts, counts, by="GENE")
    }
    index = index + 1
  }
  all_counts<-data.frame(all_counts)
  rownames(all_counts)<-all_counts$GENE
  all_counts$GENE<-NULL
  return(all_counts)
}

#read raw count dat
filelist1<-read.table(parSampleFile1, header=F, stringsAsFactors = F)
rawobjs = list()
fidx=3
fileMap<-split(filelist1$V1, filelist1$V2)

for(fileTitle in names(fileMap)) {
  fileName  = fileMap[[fileTitle]]
  cat(fileTitle, "\n")
  if(dir.exists(fileName)){
    counts = Read10X(fileName)
  } else if (grepl('.h5$', fileName)) {
    counts = Read10X_h5(fileName)
  } else if (grepl('.gz$', fileName)) {
    counts = data.frame(read_gzip_count_file(fileName, fileTitle, species))
  }
  
  adt.counts<-NULL
  if (is.list(counts) & ("Gene Expression" %in% names(counts))){
    adt.counts<-counts$`Antibody Capture`
    counts<-counts$`Gene Expression` 
  }
  
  if(fileTitle %in% names(remove_cells)){
    rcs<-remove_cells[[fileTitle]]
    if(length(rcs) > 0){
      counts<-counts[,!(colnames(counts) %in% rcs)]
      if(!is.null(adt.counts)){
        adt.counts<-adt.counts[,!(colnames(adt.counts) %in% rcs)]
      }
    }
  }

  if (species=="Mm") {
    rownames(counts)<-toMouseGeneSymbol(rownames(counts))
  }
  if (species=="Hs") {
    rownames(counts)<-toupper(rownames(counts))
  }
  sobj = CreateSeuratObject(counts = counts, project = fileTitle)
  sobj<-PercentageFeatureSet(object=sobj, pattern=Mtpattern, col.name="percent.mt")
  sobj<-PercentageFeatureSet(object=sobj, pattern=rRNApattern, col.name = "percent.ribo")
  sobj<-PercentageFeatureSet(object=sobj, pattern=hemoglobinPattern, col.name="percent.hb")    

  if (!is.null(adt.counts)){
    mat<-as.matrix(adt.counts)
    rowsum<-apply(mat>0, 1, sum)
    mat<-mat[rowsum > (ncol(mat) / 2),]
    sobj[["ADT"]] <- CreateAssayObject(counts = mat)
  }

  sobj$sample=fileTitle

  if(has_hto && fileTitle %in% names(hto_data)) {
    cell_data = hto_data[[fileTitle]]
    validobj = subset(sobj, cells=cell_data$cell)
    #tagname = unique(cell_data$HTO)[1]
    for (tagname in unique(cell_data$HTO)){
      tagcells = cell_data[cell_data$HTO == tagname,]
      sample = tagcells$Sample[1]
      sample_obj = subset(validobj, cells=tagcells$cell)
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
  if (!exists('parSampleFile3')) {
    stop('pool_sample_groups is not defined.')
  }
  if (!file.exists(parSampleFile3)) {
    stop(paste0('pool_sample_groups ', parSampleFile3, ' is not defined or not exists.'))
  }
  #cat("pooling samples ... \n")
  pools = read.table(parSampleFile3, header=F, stringsAsFactors = F)
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
rm(rawobjs)

output_rawdata(rawobj, outFile, Mtpattern, rRNApattern, hemoglobinPattern)
