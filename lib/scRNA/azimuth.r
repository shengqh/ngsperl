rm(list=ls()) 
sample_name='iSGS_scar_3364_1'
outFile='iSGS_scar_3364_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20240220_scRNA_iSGS_cell_atlas/raw_qc_sct2_Azimuth/result/iSGS_scar_3364_1')

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
azimuth_cols = c()
for(l in c("l1", "l2", "l3")){
  colname = paste0("predicted.celltype.", l)
  newcolname = paste0("Azimuth_", l)
  if(colname %in% colnames(obj@meta.data)){
    azimuth_cols = c(azimuth_cols, newcolname)
    obj <-AddMetaData(obj, metadata = obj@meta.data[,colname], col.name = newcolname)
  }else{
    colname = paste0("predicted.annotation.", l)
    if(colname %in% colnames(obj@meta.data)){
      azimuth_cols = c(azimuth_cols, newcolname)
      obj <-AddMetaData(obj, metadata = obj@meta.data[,colname], col.name = newcolname)
    }
  }
}

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))
writeLines(azimuth_cols, paste0(outFile, ".Azimuth_cols.txt"))

if("Azimuth_l2" %in% azimuth_cols){
  ct_name = "Azimuth_l2"
}else{
  ct_name = "Azimuth_l1"
}

df<-obj@meta.data[,c("orig.ident", azimuth_cols)] %>% rename(Sample=orig.ident)
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
  width=6300
}else{
  g<-g1
  width=2300
}
height=2000

ggsave(paste0(outFile, ".Azimuth.png"), g, width=width, height=height, units="px", dpi=300, bg="white")
rm(major_obj)

