rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/paula_hurley_projects/11498_WHY_VisiumHD/20251016_11498_VisiumHD_cellsegment_qc/MEcell/result/WHY_01')

### Parameter setting end ###

# Load libraries sf first, otherwise it will cause error when subset seurat object with spatial polygon assay
library(sf)

source("scRNA_func.r")
load_install("MEcell", "liuqivandy/MEcell")

myoptions_tbl=fread(parSampleFile2, header=FALSE) 
myoptions=split(myoptions_tbl$V1, myoptions_tbl$V2)

assay=myoptions$assay

is_polygons=assay == "Spatial.Polygons"
bin.size=ifelse(is_polygons, "polygons", 8)
assay_slice=ifelse(is_polygons, "slice1.polygons", "slice1.008um")

min_umi=100

data_dir <- fread(parSampleFile1, header=FALSE)$V1[1]

cat("Loading spatial data from:", data_dir, "...\n")
if(grepl("\\.rds$", tolower(data_dir))) {
  object <- readRDS(data_dir)
  DefaultAssay(object) <- assay
} else {
  object <- Seurat::Load10X_Spatial(bin.size = bin.size, data.dir = data_dir, slice = 'slice1')
}

cat(paste0("Keep the spots with at least ", min_umi, " UMIs ...\n"))
if(is_polygons){
  object <- subset(object, subset = nCount_Spatial.Polygons >= min_umi)
} else {
  object <- subset(object, subset = nCount_Spatial.008um >= min_umi)
}

cat("There are", nrow(object), "genes and", ncol(object), "spots\n")

cat("Normalizing data\n")
object <- NormalizeData(object)

cat("FindVariableFeatures\n")
object <- FindVariableFeatures(object)

cat("ScaleData\n")
object <- ScaleData(object)

cat("RunPCA\n")
object <- RunPCA(object)

cat("Run MEcell\n")
object <- MEcell(object, usepca=TRUE)

cat("FindClusters\n")
object <- FindClusters(object, graph.name="MEcell")

saveRDS(object, file=paste0(outFile, ".MEcell.rds"))

#Visualizing cell clustering
g = ImageDimPlot(object,
  group.by="MEcell_res.0.8")+ 
  ggtitle(paste0("MEcell clustering on ", assay))+
  theme(plot.title = element_text(hjust = 0.5))

if(is_polygons) {
  g = adjust_ImagePlot_Polygons(g)
}else{
  g = adjust_ImagePlot_Bins(g)
}

cluster_png=paste0(sample_name, ".MEcell.clustering.png")
ggsave(g, filename=cluster_png, width=6, height=6, dpi=300, units="in", bg="white")

