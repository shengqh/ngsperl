rm(list=ls())
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/20260120_silhouetee_pairwise/silhouette_pairwise_prepare_data/result/Aorta_Progeria.obj_data.rds'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20250513_Aorta_Progeria_scRNA_mouse/20260120_silhouetee_pairwise/silhouette_pairwise_15_16/result')

### Parameter setting end ###

library(cluster)

myoptions=read.table(parSampleFile1, header=FALSE, sep="\t")

i=as.numeric(myoptions[myoptions$V2=="i", 1])
j=as.numeric(myoptions[myoptions$V2=="j", 1])
max_cells=as.numeric(myoptions[myoptions$V2=="max_cells", 1])

cat("Reading object data file:", parFile1, "\n")

obj_data_lst=readRDS(parFile1)
clusters=obj_data_lst$clusters
coords=obj_data_lst$coords

cat("Processing cluster:", i, "and cluster:", j, "\n")
subset_indices_i <- which(clusters == i)
subset_indices_j <- which(clusters == j)

cat("There are", length(subset_indices_i), "cells in cluster", i, "and", length(subset_indices_j),"cells in cluster", j, ".\n")
if(length(subset_indices_i) > max_cells) {
  stop(paste0("Too many cells in cluster ", i, " in file ", obj_file))
}
if(length(subset_indices_j) > max_cells) {
  stop(paste0("Too many cells in cluster ", j, " in file ", obj_file))
}

get_silhouette=function(coords, clusters, cur_indices) {
  subset_coords <- coords[cur_indices, ]
  subset_labels <- clusters[cur_indices]

  cat("  Calculate distance matrix ...\n")
  dist.matrix=dist(subset_coords)

  cat("  Calculate silhouette ...\n")
  sil <- silhouette(subset_labels, dist.matrix)
  sil <- as.data.frame(sil)
  rownames(sil) <- rownames(subset_coords)

  return(sil)
}

subset_indices <- which(clusters %in% c(i, j))
if(length(subset_indices) < max_cells) {
  cat("Calculate with", length(subset_indices_i),"cluster",i,"cells and",length(subset_indices_j),"cluster",j,"cells\n")
  sil = get_silhouette(coords, clusters, subset_indices)
}else{
  cat("Because the total number of cells exceeds", max_cells, ", accurate computation of silhouette scores is computationally prohibitive; thus, we estimated silhouette scores using a subsampling-based approach.\n")
  
  set.seed(20260123)
  test_j = sample(subset_indices_j, max_cells - length(subset_indices_i))
  cat("  Calculate for cluster", i, "with", length(subset_indices_i),"cluster",i,"cells and",length(test_j),"cluster",j,"cells\n")
  cur_indices = sort(c(subset_indices_i, test_j))
  sil_i = get_silhouette(coords, clusters, cur_indices) |> dplyr::filter(cluster==i)

  set.seed(20260123)
  test_i = sample(subset_indices_i, max_cells - length(subset_indices_j))
  cat("  Calculate for cluster", j, "with", length(subset_indices_j),"cluster",j,"cells and",length(test_i),"cluster",i,"cells\n")
  cur_indices = sort(c(test_i, subset_indices_j))
  sil_j = get_silhouette(coords, clusters, cur_indices) |> dplyr::filter(cluster==j)

  sil = rbind(sil_i, sil_j)
}

sil_csv=paste0(outFile, ".", i, "_", j, ".silhouette.csv")
cat("Saving silhouette to:", sil_csv, "\n")
write.csv(sil, sil_csv, row.names=TRUE, quote=FALSE)

cat("Finished processing cluster:", i, "and cluster:", j, "\n")

