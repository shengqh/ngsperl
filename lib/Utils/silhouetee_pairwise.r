rm(list=ls()) 
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/20260120_silhouetee_pairwise/silhouette_pairwise_15_16/result')

### Parameter setting end ###

library(Seurat)
library(cluster)

myoptions=read.table(parSampleFile1, header=FALSE, sep="\t")

obj_file=myoptions[myoptions$V2==outFile, 1]
i=as.numeric(myoptions[myoptions$V2=="i", 1])
j=as.numeric(myoptions[myoptions$V2=="j", 1])
cluster=myoptions[myoptions$V2=="cluster", 1]
reduction=myoptions[myoptions$V2=="reduction", 1]

cat("Reading object file:", obj_file, "\n")
obj=readRDS(obj_file)

dims_to_use <- 1:30
coords <- Embeddings(obj, reduction = reduction)[, dims_to_use]
clusters=FetchData(obj, vars = cluster)[, 1]

cat("Processing cluster:", i, "and cluster:", j, "\n")
subset_indices <- which(clusters %in% c(i, j))
subset_coords <- coords[subset_indices, ]
subset_labels <- clusters[subset_indices]

cat("Calculate distance matrix ...\n")
dist.matrix=dist(subset_coords)

cat("Calculate silhouette ...\n")
sil <- silhouette(subset_labels, dist.matrix)
sil <- as.data.frame(sil)
rownames(sil) <- rownames(subset_coords)

sil_rds=paste0(outFile, ".", i, "_", j, ".silhouetee.rds")
cat("Saving silhouette to:", sil_rds, "\n")

result=list(i=i, j=j, silhouette=sil)
saveRDS(result, sil_rds)

cat("Finished processing cluster:", i, "and cluster:", j, "\n")

