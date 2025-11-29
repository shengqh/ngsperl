rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/paula_hurley_projects/11498_WHY_VisiumHD/20251016_11498_VisiumHD_cellsegment_qc/SpaGene/result/WHY_01')

### Parameter setting end ###

source("reportFunctions.R")
if(!requireNamespace("SpaGene", quietly = TRUE)){
  BiocManager::install("liuqivandy/SpaGene")
}

library(SpaGene)
library(Seurat)

myoptions=read_file_map(parSampleFile2, do_unlist=FALSE)
min_umi=as.numeric(myoptions$nCount_cutoff)
LRpair_file=myoptions$LRpair_file

obj_file=fread(parSampleFile1, header=FALSE)[1, V1]

cat("Read object\n")
object <- readRDS(obj_file)

# Get cell segmentation counts

DefaultAssay(object) <- "Spatial.Polygons"
counts=GetAssayData(object, assay="Spatial.Polygons", layer="counts")

location<-GetTissueCoordinates(object, image = "slice1.polygons")[,c(1,2)]

spa_rds_file=paste0(sample_name, ".SpaGene.object.rds")
if(file.exists(spa_rds_file)){
  spa<-readRDS(spa_rds_file)
}else{
  cat("SpaGene\n")
  spa<-SpaGene(counts, location)
  saveRDS(spa, file=spa_rds_file)
}

spa_pattern_rds=paste0(sample_name, ".SpaGene.pattern.rds")
if(file.exists(spa_pattern_rds)){
  spa_pattern<-readRDS(spa_pattern_rds)
}else{
  cat("SpaGene_Pattern\n")
  spa_pattern<-FindPattern(spa,nPattern = 10)
  saveRDS(spa_pattern, file=spa_pattern_rds)
}

spa_lr_rds=paste0(sample_name, ".SpaGene.LR.rds")
if(file.exists(spa_lr_rds)){
  spa_lr<-readRDS(spa_lr_rds)
}else{
  cat("SpaGene_LR\n")
  load(LRpair_file)
  spa_lr<-SpaGene_LR(counts,location,LRpair=LRpair)  
  saveRDS(spa_lr, file=spa_lr_rds)
}

cat("Done.\n")
