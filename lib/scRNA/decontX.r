rm(list=ls()) 
sample_name='Adipose_9240'
outFile='Adipose_9240'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20241224_combined_scRNA_hg38_cellbender_fastmnn/cellbender_raw_qc_decontX/result/Adipose_9240')

### Parameter setting end ###

source("scRNA_func.r")
load_install("SingleCellExperiment")
load_install("celda")
load_install("ggplot2")
load_install("patchwork")

sample_map = read_file_map("fileList1.txt")
raw_map = read_file_map("fileList2.txt")

myoptions = read_file_map("fileList3.txt", do_unlist=FALSE)
myoptions$remove_decontX = is_one(myoptions$remove_decontX)
myoptions$remove_decontX_by_contamination = as.numeric(myoptions$remove_decontX_by_contamination)

read_sce<-function(countfile, species, is_qc_object=FALSE, sample_name=""){
  cat("  read", countfile, "\n")
  if(is_qc_object){
    obj<-readRDS(countfile)
    if(sample_name %in% names(obj)){
      obj=obj[[sample_name]]
    }
    if(is.list(obj)){
      obj=obj$obj
    }
    counts = MyGetAssayData(obj, "RNA", "counts")
  }else{
    lst = read_scrna_data(countfile)
    counts<-lst$counts
  }

  if (species=="Mm") {
    rownames(counts)<-toMouseGeneSymbol(rownames(counts))
  }
  if (species=="Hs") {
    rownames(counts)<-toupper(rownames(counts))
  }

  sce <- SingleCellExperiment(assays=list(counts=counts))
  return(sce)
}

draw_umap<-function(sce, png_file){
  df = data.frame(sce@metadata$decontX$estimates$all_cells$UMAP)
  df$contamination = unlist(sce@metadata$decontX$estimates$all_cells$contamination)

  #if we subset sce, metadata would not be changed, so we need to subset df as well
  df = df[colnames(sce),]
  df = df[order(df$contamination),]

  g1<-ggplot(df, aes(DecontX_UMAP_1, DecontX_UMAP_2)) + geom_point(aes(color=contamination), size=0.1) + scale_color_gradient(low="lightgray", high="red", limits = c(0,1)) + theme_bw3() + theme(aspect.ratio = 1)
  g2<-ggplot(df, aes(x = contamination)) + geom_histogram(aes(y = after_stat(density)), color="black", fill=NA, bins = 50) + geom_density() + theme_bw3() + xlim(0,1)
  g<-g1+g2+plot_layout(ncol=2)
  png(png_file, width=2300, height=1000, res=300)
  print(g)
  dev.off()
}

sample_name=outFile
cat(sample_name, "\n")
countfile = sample_map[sample_name]
raw_countfile = raw_map[sample_name]

clusters = NULL
is_qc_object=FALSE
if(grepl(".rds$", countfile)){
  #countfile is a object file, containing cell_type/cluster information
  obj<-readRDS(countfile)
  if(sample_name %in% names(obj)){
    obj=obj[[sample_name]]
    is_qc_object=TRUE
  }
  if(is.list(obj)){
    obj=obj$obj
  }
  if("layer4" %in% colnames(obj@meta.data)){
    clusters = obj@meta.data$layer4
  }else if("cell_type" %in% colnames(obj@meta.data)){
    clusters = obj@meta.data$cell_type
  }
}

sce <- read_sce(countfile, myoptions$species, is_qc_object=is_qc_object, sample_name=sample_name)
sce.raw <- read_sce(raw_countfile, myoptions$species, is_qc_object=FALSE, sample_name="")

common_genes = intersect(row.names(sce), row.names(sce.raw))
sce=sce[common_genes,]
sce.raw=sce.raw[common_genes,]

sce <- decontX(sce, background = sce.raw, z=clusters)

draw_umap(sce, paste0(sample_name, ".decontX.png"))

meta=colData(sce)
saveRDS(meta, paste0(sample_name, ".decontX.meta.rds"))

if(myoptions$remove_decontX & (myoptions$remove_decontX_by_contamination > 0)){
  cat("  remove cells with contamination > ", myoptions$remove_decontX_by_contamination, "\n")
  decontX_sce = sce[, sce@metadata$decontX$estimates$all_cells$contamination < myoptions$remove_decontX_by_contamination]

  draw_umap(decontX_sce, paste0(sample_name, ".decontX.after.png"))

  filtered = data.frame("Pre_filter_cell" = ncol(sce), "Post_filter_cell" = ncol(decontX_sce))
  write.csv(filtered, paste0(sample_name, ".decontX.filtered.csv"), row.names=FALSE)

  #we only remove the cells but not change the gene read counts
  counts = counts(decontX_sce)
}else{
  #update the gene read count
  counts = ceiling(decontXcounts(sce))
}
saveRDS(counts, paste0(sample_name, ".decontX.counts.rds"))

