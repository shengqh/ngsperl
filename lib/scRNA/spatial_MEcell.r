rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20250620_spatial_scRNA_analysis/MEcell/result/WHY_01')

### Parameter setting end ###

source("scRNA_func.r")
library(MEcell)

myoptions_tbl=fread(parSampleFile2, header=FALSE) 
myoptions=split(myoptions_tbl$V1, myoptions_tbl$V2)

assay=myoptions$assay
sketch_assay=myoptions$sketch_assay

obj_file=fread(parSampleFile1, header=FALSE) |>
  dplyr::filter(V2 == sample_name) |>
  dplyr::pull(V1)

cat("Read object from", obj_file, "...\n")

object <- readRDS(obj_file)
DefaultAssay(object) <- assay
cat("There are ", nrow(object), "genes and", ncol(object), "bins in the object\n")

cat("Keep sketched bins only ...")
sketched_counts=GetAssayData(object, assay=sketch_assay, layer="counts")
object=subset(object, cells=colnames(sketched_counts))
object@assays[[sketch_assay]] = NULL

cat("There are ", nrow(object), "genes and", ncol(object), "bins in the sketched object\n")

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

#Visualizing cell clustering
g = ImageDimPlot(object,
  group.by="MEcell_res.0.8")+ 
  ggtitle("MEcell clustering")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_reverse()

saveRDS(object, file=paste0(outFile, ".MEcell.rds"))

cluster_png=paste0(sample_name, ".MEcell.clustering.png")
ggsave(g, filename=cluster_png, width=6, height=6, dpi=300, units="in", bg="white")

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
