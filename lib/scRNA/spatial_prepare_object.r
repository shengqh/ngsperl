rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20250414_T01_prepare_bin_counts/prepare_spatial_object/result/WHY_01')

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
bin.size=8
DefaultAssay = "Spatial.008um"
Cutoff<-data.frame(nFeature_cutoff_min = 20 ,nFeature_cutoff_max = 1000,nCount_cutoff=40, mt_cutoff = 25, cluster_remove=c(""),stringsAsFactors = F)

###########################
# message
###########################
print(paste0("Working in ",data_path))
print(paste0("sample_name: ",sample_name,"; bin.size: ",bin.size,"; DefaultAssay: ",DefaultAssay))

# 1.Load and processing Objects _________________________________________
###########################
# Load Objects
###########################
object <- Load10X_Spatial(data.dir = data_path, bin.size = c(bin.size))
object$sample <- sample_name

###########################
# Default Assay 8µm
###########################
Assays(object)
DefaultAssay(object) <- DefaultAssay

###########################
# Mitochondrial Percentage
###########################
object[['percent.mt']] <- PercentageFeatureSet(object, pattern = "^MT-")
object <- subset(object, subset = nFeature_Spatial.008um > Cutoff$nFeature_cutoff_min & nFeature_Spatial.008um<Cutoff$nFeature_cutoff_max & nCount_Spatial.008um>Cutoff$nCount_cutoff & percent.mt < Cutoff$mt_cutoff)

###########################
# Normalize
###########################
object <- NormalizeData(object)

# 3. CLUSTERING___________________________________________
###########################
# Variable Feature and scaling
###########################
object <- FindVariableFeatures(object)
object <- ScaleData(object)

# we select 5,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "RNA"
)

# 4. Save rds object___________________________________________
###########################
# save rds
###########################
saveRDS(object, file = paste0(sample_name, "_",bin.size,"um.rds"))
