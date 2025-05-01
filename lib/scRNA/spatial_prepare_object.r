rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20250415_spatial_scRNA_analysis/individual_spatial_object/result/WHY_01')

### Parameter setting end ###

source("scRNA_func.r")
set.seed(20250414)

library(Seurat)

plan("multicore", workers = 16)
options(future.globals.maxSize = 10000 * 1024^2)
options(future.seed = TRUE)

file_map=read_file_map(parSampleFile1, do_unlist = FALSE)
data_path = file_map[[sample_name]]

###########################
# Parameters
###########################
myoptions_tbl=fread('fileList2.txt', header = F, sep = '\t')
myoptions=split(myoptions_tbl$V1, myoptions_tbl$V2)
bin.size=as.numeric(myoptions$bin.size)
if(bin.size == 8){
  binstr="008um"
}else{
  stop("Not implemented for bin size other than 8um")
}

assay = paste0("Spatial.", binstr)

nFeature_cutoff_min = as.numeric(myoptions$nFeature_cutoff_min)
nFeature_cutoff_max = as.numeric(myoptions$nFeature_cutoff_max)
nCount_cutoff = as.numeric(myoptions$nCount_cutoff)
mt_cutoff = as.numeric(myoptions$mt_cutoff)

Cutoff<-data.frame( nFeature_cutoff_min = nFeature_cutoff_min ,
                    nFeature_cutoff_max = nFeature_cutoff_max,
                    nCount_cutoff=nCount_cutoff, 
                    mt_cutoff = mt_cutoff, 
                    cluster_remove=c(""),
                    stringsAsFactors = F)

cat("Working in",data_path, "\n")
cat("sample_name:",sample_name,"; bin.size:",bin.size,"; DefaultAssay:",assay, "\n")

cat("Loading data...\n")
obj <- Load10X_Spatial(data.dir = data_path, bin.size = c(bin.size))
obj$orig.ident <- sample_name
obj$sample <- sample_name

Assays(obj)
DefaultAssay(obj) <- assay

obj<-PercentageFeatureSet(object=obj, pattern=myoptions$Mtpattern, col.name = "percent.mt")
obj<-PercentageFeatureSet(object=obj, pattern=myoptions$rRNApattern, col.name = "percent.ribo")

cells=colnames(obj)[obj@meta.data[, paste0("nFeature_", assay)] > nFeature_cutoff_min & 
                    obj@meta.data[, paste0("nFeature_", assay)] < nFeature_cutoff_max &
                    obj@meta.data[, paste0("nCount_", assay)] > nCount_cutoff &
                    obj@meta.data[, "percent.mt"] < mt_cutoff]
cat("Cells before filter=", ncol(obj), "\n")
obj <- subset(obj, cells = cells)
cat("Cells after filter=", ncol(obj), "\n")

cat("Normalizing data...\n")
obj <- NormalizeData(obj)

cat("Finding variable features...\n")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

cat("Save object...\n")
saveRDS(obj, file = paste0(sample_name, "_", bin.size, "um.full.rds"))

cat("SketchData...\n")
#https://github.com/satijalab/seurat/issues/7171
#without features = VariableFeatures(obj) , it will use all features and throw error: too slow.
obj <- SketchData(
  obj = obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "RNA",
  features = VariableFeatures(obj) 
)

obj[[assay]] <- NULL
counts=GetAssayData(obj, assay = "RNA", layer = "counts")
obj@meta.data=obj@meta.data[colnames(counts),,drop=FALSE]
obj@commands=list()

cat("Save sketch object...\n")
saveRDS(obj, file = paste0(sample_name, "_", bin.size, "um.sketch.rds"))

cat("Done.\n")

