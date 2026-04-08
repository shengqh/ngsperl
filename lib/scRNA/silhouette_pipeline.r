rm(list=ls()) 
outFile='Aorta_Progeria'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_sct2_fastmnn_dr0.5_3_choose/result/Aorta_Progeria.final.rds'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_sct2_fastmnn_dr0.5_3_choose_silhouette_pca/result')

### Parameter setting end ###

obj_file <- parFile1

myoptions_tbl=read.table(parSampleFile1, header=FALSE, stringsAsFactors=FALSE)
myoptions=split(myoptions_tbl$V1, myoptions_tbl$V2)

silhouette_bin <- myoptions$silhouette_path
cluster_col <- myoptions$cluster_col
reduction <- myoptions$reduction

cat("obj_file:", obj_file, "\n")
cat("cluster_col:", cluster_col, "\n")
cat("reduction:", reduction, "\n")
cat("outFile:", outFile, "\n")
cat("silhouette_bin:", silhouette_bin, "\n")

library(Seurat)
library(data.table)

seurat_coords_file <- paste0(outFile, ".seurat.coords.csv")
seurat_clusters_file <- paste0(outFile, ".seurat.clusters.txt")
seurat_sil_file <- paste0(outFile, ".seurat.silhouette.csv")
consistent_clusters_file <- paste0(outFile, ".consistent.clusters.txt")
consistent_sil_file <- paste0(outFile, ".consistent.silhouette.csv")
closest_file <- paste0(outFile, ".silhouette.closest.csv")

# Step 1: Extract coordinates and clusters from seurat object
cat("Step 1: Extracting coordinates and clusters from seurat object\n")
cat("Reading object file:", obj_file, "\n")
obj <- readRDS(obj_file)
clusters <- FetchData(obj, vars = cluster_col)[, 1]
coords <- Embeddings(obj, reduction = reduction)
writeLines(as.character(clusters), seurat_clusters_file)
write.csv(coords, seurat_coords_file, row.names=TRUE, quote=FALSE)
cat("Saved clusters to:", seurat_clusters_file, "\n")
cat("Saved coords to:", seurat_coords_file, "\n")
rm(obj)

# Step 2: Calculate silhouette for seurat clusters
cat("Step 2: Calculating silhouette for seurat clusters\n")
cmd <- paste(shQuote(silhouette_bin), shQuote(seurat_coords_file), shQuote(seurat_clusters_file), shQuote(seurat_sil_file))
cat("Running:", cmd, "\n")
ret <- system(cmd)
if (ret != 0) stop("silhouette command failed with exit code ", ret)

# Step 3: Set inconsistent cells to cluster -1
cat("Step 3: Setting inconsistent cells to cluster -1\n")
seurat_sil_data <- fread(seurat_sil_file)
coords_data <- data.frame(fread(seurat_coords_file), row.names=1)
stopifnot(all(rownames(coords_data) %in% seurat_sil_data$name))
seurat_sil_data$new_cluster <- ifelse(seurat_sil_data$sil_width > 0, seurat_sil_data$cluster, -1)
writeLines(as.character(seurat_sil_data$new_cluster), consistent_clusters_file)
cat("Saved consistent clusters to:", consistent_clusters_file, "\n")

# Step 4: Calculate silhouette for consistent clusters
cat("Step 4: Calculating silhouette for consistent clusters\n")
cmd <- paste(shQuote(silhouette_bin), shQuote(seurat_coords_file), shQuote(consistent_clusters_file), shQuote(consistent_sil_file))
cat("Running:", cmd, "\n")
ret <- system(cmd)
if (ret != 0) stop("silhouette command failed with exit code ", ret)

# Step 5: Combine results - assign inconsistent cells to closest cluster
cat("Step 5: Combining silhouette results\n")
seurat_sil_data <- fread(seurat_sil_file)
consistent_sil_data <- fread(consistent_sil_file)

stopifnot(all(consistent_sil_data$name == seurat_sil_data$name))

consistent_sil_data$seurat_cluster <- seurat_sil_data$cluster

consistent_sil_data$closest <- ifelse(consistent_sil_data$sil_width > 0,
                                      consistent_sil_data$seurat_cluster,
                                      consistent_sil_data$neighbor)

result <- consistent_sil_data |>
  dplyr::select(name, seurat_cluster, closest, sil_width) |>
  dplyr::rename(cluster = seurat_cluster)

write.csv(result, closest_file, row.names=FALSE, quote=FALSE)
cat("Saved closest clusters to:", closest_file, "\n")

