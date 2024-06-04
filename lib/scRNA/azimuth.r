rm(list=ls()) 
sample_name='iSGS_scar_3364_1'
outFile='iSGS_scar_3364_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240320_scRNA_iSGS_cell_atlas/raw_qc_sct2_Azimuth/result/iSGS_scar_3364_1')

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

Amimuth_ref = myoptions$Azimuth_ref

if(!exists("obj")){
  obj=read_object_from_file_list(parSampleFile1)
}

obj <- RunAzimuth(query = obj, reference = Amimuth_ref)

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
