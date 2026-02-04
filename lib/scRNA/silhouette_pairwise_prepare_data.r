rm(list=ls()) 
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/20260120_silhouetee_pairwise/silhouette_pairwise_01_prepare_data/result')

### Parameter setting end ###

library(Seurat)
library(cluster)

myoptions=read.table(parSampleFile1, header=FALSE, sep="\t")

obj_file=myoptions[myoptions$V2==outFile, 1]
cluster=myoptions[myoptions$V2=="cluster", 1]
reduction=myoptions[myoptions$V2=="reduction", 1]

cat("Reading object file:", obj_file, "\n")
obj=readRDS(obj_file)

clusters=FetchData(obj, vars = cluster)[, 1]

dims_to_use <- 1:30
coords <- Embeddings(obj, reduction = reduction)[, dims_to_use]

obj_data=list(clusters=clusters, coords=coords)

data_rds=paste0(outFile, ".obj_data.rds")
cat("Saving clusters and coords to:", data_rds, "\n")
saveRDS(obj_data, data_rds)

