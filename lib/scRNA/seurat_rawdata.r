
source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(digest)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
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
    cell_data = read.csv(cell_file, stringsAsFactors=FALSE, header=TRUE, check.names=F)
    colnames(cell_data)[1]<-"cell"
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
    sample_obj = object.list[[sample_id]]
    cd<-FetchData(sample_obj, c("seurat_clusters"))
    rm(sample_obj)
    
    crs<-objs_config$cluster_remove[rowi]
    crslist<-unlist(strsplit(crs,','))
    rcs<-rownames(cd)[as.character(cd$seurat_clusters) %in% crslist] 
    
    remove_cells[[sample_id]]=rcs
  }
}

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
  sobj[["percent.mt"]] <- PercentageFeatureSet(object = sobj, pattern = Mtpattern)
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

writeLines(rownames(rawobj), paste0(outFile, ".genes.txt"))

saveRDS(rawobj, paste0(outFile, ".rawobj.rds"));

png(file=paste0(outFile, ".qc.png"), width=3000, height=1200, res=300)
p1 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(object = rawobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p<-p1+p2+plot_layout(ncol=2)
print(p)
dev.off()

nsample=length(unique(rawobj$sample))
mt<-data.frame(mt=rawobj$percent.mt, Sample=rawobj$sample, nFeature=log10(rawobj$nFeature_RNA), nCount=log10(rawobj$nCount_RNA))
png(file=paste0(outFile, ".qc.individual.png"), width=3000, height=min(20000, 1200 * nsample), res=300)
p1<-ggplot(mt, aes(x=mt,y=nCount) ) +
  geom_bin2d(bins = 70) + 
  scale_fill_continuous(type = "viridis") + 
  xlab("Percentage of mitochondrial") + ylab("log10(number of read)") +
  facet_grid(Sample~.) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
p2<-ggplot(mt, aes(x=mt,y=nFeature) ) +
  geom_bin2d(bins = 70) + 
  scale_fill_continuous(type = "viridis") + 
  xlab("Percentage of mitochondrial") + ylab("log10(number of feature)") +
  facet_grid(Sample~.) + theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"))
p<-p1+p2+plot_layout(ncol=2)
print(p)
dev.off()
