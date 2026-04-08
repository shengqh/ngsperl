rm(list=ls())
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_fastmnn_dr0.5_3_choose/result/Aorta_Progeria.final.rds'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_fastmnn_dr0.5_3_choose_silhouette_1_prepare_data/result')


### Parameter setting end ###

library(Seurat)
library(cluster)

myoptions=read.table(parSampleFile1, header=FALSE, sep="\t")

obj_file=parFile1
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

writeLines(as.character(clusters), paste0(outFile, ".clusters.txt"))
write.csv(coords, paste0(outFile, ".coords.csv"), row.names=TRUE, quote=FALSE)
