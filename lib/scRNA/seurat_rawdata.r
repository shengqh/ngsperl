rm(list=ls()) 
outFile='GPA'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20230807_gpa_scRNA_hg38/decontX_nd_seurat_rawdata.test/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(digest)
library(patchwork)
library(sparseMatrixStats)
library(data.table)
library(tidyr)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

myoptions<-read_file_map(parSampleFile2, do_unlist = FALSE)

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
hemoglobinPattern=myoptions$hemoglobinPattern

species=myoptions$species
pool_sample<-is_one(myoptions$pool_sample)

keep_seurat_object<-is_one(myoptions$keep_seurat_object)
seurat_sample_column<-myoptions$seurat_sample_column

ensembl_map=NULL
if("ensembl_gene_map_file" %in% names(myoptions)){
  ensembl_gene_map_file = myoptions$ensembl_gene_map_file
  if(ensembl_gene_map_file != ""){
    gene_tb=read.table(ensembl_gene_map_file, sep="\t", header=T)
    gene_tb=gene_tb[!duplicated(gene_tb$ENSEMBL_GENE_ID),]
    gene_tb=gene_tb[gene_tb$ENSEMBL_GENE_ID != "",]
    ensembl_map = split(gene_tb$GENE_SYMBOL, gene_tb$ENSEMBL_GENE_ID)
  }
}

hto_str = ""
hto_data = list()
hto_sample_file<-myoptions$hto_sample_file
has_hto = (hto_sample_file != "")

if (has_hto) {
  if ("hto_ignore_samples" %in% names(myoptions)){
    hto_ignore_samples = myoptions$hto_ignore_samples
  }else{
    hto_ignore_samples = ""
  }
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

read_remove_cells_from_combined_qc<-function(obj_file, config_file){
  result = list()

  load(obj_file)
  
  objs_config=read.table(config_file, sep="\t", header=T, stringsAsFactors = F)
  
  if(all(is.null(names(object.list)))){
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
        
        result[[sample_id]]=rcs
      }
    }
  }
  return(result)
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

#read object from either qc seurat object file or original count file
raw_objs = list()
remove_cells = list()

is_qc_data = exists('parSampleFile6')
if(is_qc_data){
  #from qc data
  obj_map = read_file_map(parSampleFile6)
  if(parFile1 != ""){
    remove_tbl = fread(parFile1, header=T, data.table=F)
    remove_tbl = remove_tbl[remove_tbl$cluster_remove != "",,drop=F]
    if(nrow(remove_tbl) == 0){
      remove_map = list()
    }else{
      remove_tbl = separate_rows(remove_tbl, cluster_remove, sep=';')
      remove_map = split(remove_tbl$cluster_remove, remove_tbl$sample)
    }
  }else{
    remove_map = list()
  }

  sample_names = names(obj_map)
  sample_name = sample_names[1]
  for(sample_name in sample_names){
    obj_file = obj_map[[sample_name]]
    cat("reading", sample_name, ":", obj_file, "\n")
    sobj = read_object(obj_file, sample_name=sample_name)

    if(sample_name %in% names(remove_map)){
      cd<-sobj@meta.data
      crslist<-as.character(remove_map[[sample_name]])
      rm_cells = rownames(cd)[as.character(cd$seurat_clusters) %in% crslist] 
      remove_cells[[sample_name]] = rm_cells
      cat("  ", length(rm_cells), "cells from clusters", paste0(crslist, collapse=", "), "will be removed.\n")
    }
    
    raw_objs[[sample_name]]=sobj
  }
}else{
  if(file.exists(parFile1) & file.exists(parFile2)){
    remove_cells = read_remove_cells_from_combined_qc(parFile1, parFile2)
  }

  if(exists("parSampleFile5")){
    doublet_tbl=read.table(parSampleFile5, sep="\t", header=F)
    doublet_map=unlist(split(doublet_tbl$V1, doublet_tbl$V2))
    doublet_column=myoptions$doublet_column
    doublet_cells = list()
    sname = names(doublet_map)[1]
    for(sname in names(doublet_map)){
      sfile = doublet_map[sname]
      smeta = readRDS(sfile)
      dcells = rownames(smeta)[smeta[,doublet_column] %in% c("Doublet", "doublet")]
      if(sname %in% names(remove_cells)){
        remove_cells[[sname]] = c(remove_cells[[sname]], dcells)
      }else{
        remove_cells[[sname]] = dcells
      }
    }
    remove_doublets=TRUE
  }


  #read raw count dat
  file_map<-read_file_map(parSampleFile1)
  if(length(file_map) == 0){
    die("No data defined, check your configuration file.")
  }

  sample_names = names(file_map)
  sample_name = sample_names[1]
  for(sample_name in sample_names) {
    file_path  = file_map[[sample_name]]
    cat("reading", sample_name, ":", file_path, "\n")
    lst = read_scrna_data(file_path, keep_seurat_object)

    counts = lst$counts  
    adt.counts = lst$adt.counts

    if(is_seurat_object(counts)){
      sobj<-counts
      rm(counts)

      if("" != seurat_sample_column){
        sobj$orig.ident = unlist(sobj@meta.data[,seurat_sample_column])
      }else{
        sobj$orig.ident = sample_name
      }
    }else{
      rs<-rowSums(counts)
      counts<-counts[rs>0,]

      if(!is.null(ensembl_map)){
        gtf_counts<-counts[!(rownames(counts) %in% names(ensembl_map)),]
        ensembl_counts<-counts[(rownames(counts) %in% names(ensembl_map)),]
        gene_names<-unlist(ensembl_map[rownames(ensembl_counts)])
        gene_counts<-DelayedArray::rowsum(ensembl_counts, gene_names)
        counts<-rbind(gtf_counts, gene_counts)
      }

      if (species=="Mm") {
        rownames(counts)<-toMouseGeneSymbol(rownames(counts))
      }
      if (species=="Hs") {
        rownames(counts)<-toupper(rownames(counts))
      }
      sobj = CreateSeuratObject(counts = counts, project = sample_name)
      sobj$orig.ident <- sample_name
      rm(counts)
    }

    if (!is.null(adt.counts)){
      mat<-as.matrix(adt.counts)
      rowsum<-apply(mat>0, 1, sum)
      mat<-mat[rowsum > (ncol(mat) / 2),,drop=FALSE]
      if(nrow(mat) > 0){
        sobj[["ADT"]] <- CreateAssayObject(counts = mat)
      }
    }

    raw_objs[[sample_name]]=sobj
  }
}

sample_names=names(raw_objs)
sample_name=sample_names[1]
for(sample_name in sample_names){
  sobj = raw_objs[[sample_name]]

  if(!("orig.cell" %in% colnames(sobj@meta.data))){
    sobj$orig.cell = colnames(sobj)
  }

  if(!is_qc_data){
    sobj<-PercentageFeatureSet(object=sobj, pattern=Mtpattern, col.name="percent.mt", assay="RNA")
    sobj<-PercentageFeatureSet(object=sobj, pattern=rRNApattern, col.name = "percent.ribo", assay="RNA")
    sobj<-PercentageFeatureSet(object=sobj, pattern=hemoglobinPattern, col.name="percent.hb", assay="RNA")
  }

  if(sample_name %in% names(remove_cells)){
    rcs<-remove_cells[[sample_name]]
    if(length(rcs) > 0){
      valid_cells=colnames(sobj)[!(colnames(sobj) %in% rcs)]
      sobj=subset(sobj, cells=valid_cells)
    }
  }

  #################################################
  #three columns for source of the data
  #project indicates the h5 file name (experimental name)
  #sample indicates the sample name (after demultiplex)
  #orig.ident indicates final name (after pooling)
  #################################################
  #  with demultiplex: 
  #    sample != project
  #    with pooling:
  #      orig.ident != sample
  #    without pooling:
  #      orig.ident == sample
  #  without demultiplex:
  #    sample == project
  #    with pooling:
  #      orig.ident != sample
  #    without pooling:
  #      orig.ident == sample
  #################################################
  sobj$project=sample_name    
  sobj$sample=sample_name

  b_rename_cells = all(!startsWith(colnames(sobj), sample_name))
  if(sample_name %in% names(hto_data)) {
    cat("processing HTO demultiplex of", sample_name, "\n")
    raw_objs[[sample_name]] = NULL
    cell_data = hto_data[[sample_name]]
    validobj = subset(sobj, cells=cell_data$cell)
    #tagname = unique(cell_data$HTO)[1]
    for (tagname in unique(cell_data$HTO)){
      tagcells = cell_data[cell_data$HTO == tagname,]
      sample = tagcells$Sample[1]
      if(sample %in% hto_ignore_samples){
        cat("skipping HTO sample", sample, "in file", sample_name, "\n")
        next
      }

      sample_obj = subset(validobj, cells=tagcells$cell)
      sample_obj$orig.ident = sample
      sample_obj$sample = sample

      if(b_rename_cells){
        sample_obj<-RenameCells(object=sample_obj, new.names=paste0(sample_name, "_", colnames(sample_obj)))
      }

      raw_objs[[sample]] = sample_obj
    }
  }else{
    if(b_rename_cells){
      sobj<-RenameCells(object=sobj, new.names=paste0(sample_name, "_", colnames(sobj)))
    }
    raw_objs[[sample_name]] = sobj
  }
}

if(pool_sample){
  if (!exists('parSampleFile3')) {
    stop('pool_sample_groups is not defined.')
  }
  if (!file.exists(parSampleFile3)) {
    stop(paste0('pool_sample_groups ', parSampleFile3, ' is not defined or not exists.'))
  }
  cat("pooling samples ... \n")
  pool_tbl = read.table(parSampleFile3, header=F, stringsAsFactors = F)
  pool_map = split(pool_tbl$V1, pool_tbl$V2)
  pool_names = names(pool_map)
  pool_name = pool_names[1]
  for(pool_name in pool_names){
    s_names = pool_map[[pool_name]]
     #Checking for agreement in names
    if(!all(s_names %in% names(raw_objs))){
      missed = s_names[!(s_names %in% names(raw_objs))]
      stop(paste0(paste0(missed, collapse=","), " of ", pool_name, " not in object list: ", paste0(names(raw_objs), collapse=",")))
    }
    
    curListIndecies<-which(names(raw_objs) %in% s_names)

    cur_objs<-raw_objs[curListIndecies]
    cur_names = names(raw_objs)[curListIndecies]
    if(length(cur_objs) == 1){
      cur_obj = cur_objs[[1]]
    }else{
      cur_obj <- merge(cur_objs[[1]], y = unlist(cur_objs[2:length(cur_objs)]), project = "integrated")
    }
    cur_obj$orig.ident=pool_name
    Idents(cur_obj)<-"orig.ident"

    for(cur_name in cur_names){
      raw_objs[[cur_name]] = NULL
    }
    raw_objs[[pool_name]] = cur_obj
  }
}

if(length(raw_objs) == 1){
  rawobj <- raw_objs[[1]]
}else{
  cat("merging all samples ... \n")
  rawobj <- merge(raw_objs[[1]], y = unlist(raw_objs[2:length(raw_objs)]), project = "integrated")
}
rm(raw_objs)

cat("outputing result ... \n")

saveRDS(rawobj, paste0(outFile, ".rawobj.rds"))

#rawobj <-readRDS(paste0(outFile, ".rawobj.rds"))

output_rawdata(rawobj, outFile, Mtpattern, rRNApattern, hemoglobinPattern)
