rm(list=ls()) 
sample_name='KrasnowHLCA_P3_8'
outFile='KrasnowHLCA_P3_8'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20241217_endothelial_isgs_lung/20241217_endothelial_lung/raw_qc_Azimuth/result/KrasnowHLCA_P3_8')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

if(file.exists(parSampleFile3)) {
  Azimuth_ref_tbl = fread(parSampleFile3, header=F, stringsAsFactors = F)
  Azimuth_ref_map = split(Azimuth_ref_tbl$V1, Azimuth_ref_tbl$V2)
  Azimuth_ref = Azimuth_ref_map[[sample_name]]
}else{
  Azimuth_ref=NULL
}

if(is.null(Azimuth_ref)){
  Azimuth_ref=myoptions$Azimuth_ref
}

if(is.null(Azimuth_ref)){
  stop("Azimuth_ref is not defined")
}

cat("Azimuth_ref: ", Azimuth_ref, "\n")

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

obj <- RunAzimuth(query = obj, reference = Azimuth_ref)

anno_columns=grep('predicted.+\\d$', colnames(obj@meta.data), value=TRUE)
azimuth_cols = c()
idx=1
for(idx in 1:length(anno_columns)){
  colname = anno_columns[idx]
  newcolname = paste0("Azimuth_l", idx)
  azimuth_cols = c(azimuth_cols, newcolname)
  obj <-AddMetaData(obj, metadata = obj@meta.data[,colname], col.name = newcolname)
}
writeLines(azimuth_cols, paste0(outFile, ".Azimuth_cols.txt"))

if(any(grepl("predicted.+finest", colnames(obj@meta.data)))) {
  anno_columns=grep("predicted.+finest", colnames(obj@meta.data), value=TRUE)
  finest_column=unique(gsub('.score$','',anno_columns))[1]
}else{
  if("Azimuth_l2" %in% azimuth_cols){
    finest_column = "Azimuth_l2"
  }else{
    finest_column = "Azimuth_l1"
  }
}
ct_name = "Azimuth_finest"
obj <-AddMetaData(obj, metadata=obj@meta.data[,finest_column], col.name=ct_name)

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

df<-obj@meta.data[,c("orig.ident", ct_name)] %>% rename(Sample=orig.ident)
df_tbl<-table(df[,ct_name],df$Sample)
write.csv(df_tbl, paste0(outFile, ".Azimuth_Sample.csv"))

major_obj<-subset(obj, cells=colnames(obj)[!is.na(obj@meta.data[,ct_name])])
rm(obj)

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

g1=MyDimPlot(major_obj, group.by = ct_name, reduction="umap", label=T)

if(has_bubblemap){
  g2<-get_bubble_plot(major_obj, NA, ct_name, bubblemap_file, assay="RNA", 
    species=myoptions$species)
  layout <- "ABB"
  g<-g1+g2+plot_layout(design=layout)
  width=10000
}else{
  g<-g1
  width=4300
}
height=2000

ggsave(paste0(outFile, ".Azimuth.png"), g, width=width, height=height, units="px", dpi=300, bg="white")
rm(major_obj)

if(dir.exists(".local")){
  unlink(".local", recursive=TRUE)
}
