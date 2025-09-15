rm(list=ls()) 
sample_name='S12_SolidPTC_NCOA4RET'
outFile='S12_SolidPTC_NCOA4RET'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20250812_VisiumHD_RCTD_test/RCTD_singlet/result/S12_SolidPTC_NCOA4RET')

### Parameter setting end ###

library(Seurat)
library(dplyr)
library(data.table)

options_df=fread("fileList2.txt", header=FALSE)
myoptions=split(options_df$V1, options_df$V2)

obj_file <- fread(parSampleFile1, header=FALSE)$V1[1]
obj <- readRDS(obj_file)

singlet_obj=subset(obj, subset = RCTD1_Class == "singlet")
singlet_counts <- GetAssayData(singlet_obj, assay="Spatial.008um", layer="counts")

rna_obj=CreateSeuratObject(singlet_counts)
rna_obj@meta.data$orig.ident=sample_name
rna_obj <- AddMetaData(rna_obj, metadata = singlet_obj@meta.data |> dplyr::select(-orig.ident, -nCount_Spatial.008um, -nFeature_Spatial.008um))

saveRDS(rna_obj, paste0(sample_name, ".post_RCTD.singlet.RDS"))
