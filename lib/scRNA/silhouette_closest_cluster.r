args <- commandArgs(trailingOnly = TRUE)

if(0){
  args=c( "/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_sct2_fastmnn_dr0.5_3_choose_silhouette_fastmnn/result/Aorta_Progeria.fastmnn.seurat.silhouette.csv", 
          "/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_sct2_fastmnn_dr0.5_3_choose_silhouette_fastmnn/result/Aorta_Progeria.fastmnn.consistent.silhouette.csv",
          "/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_sct2_fastmnn_dr0.5_3_choose_silhouette_fastmnn/result/Aorta_Progeria.fastmnn.silhouette.closest.csv")
}

if (length(args) != 3) {
  stop("Usage: Rscript silhouette_data_from_sihouette.r <silhouette_file> <coords_file> <output_prefix>")
}

seurat_silhouette_file <- args[1]
consistent_silhouette_file <- args[2]
outFile <- args[3]

cat("seurat_silhouette_file:", seurat_silhouette_file, "\n")
cat("consistent_silhouette_file:", consistent_silhouette_file, "\n")
cat("outFile:", outFile, "\n")

library(data.table)

cat("Reading seurat silhouette file:", seurat_silhouette_file, "\n")
seurat_silhouette_data <- fread(seurat_silhouette_file) |>
  dplyr::select(name, cluster) |>
  dplyr::rename(seurat_cluster = cluster)

cat("Reading consistent silhouette file:", consistent_silhouette_file, "\n")
consistent_silhouette_data <- fread(consistent_silhouette_file)

stopifnot(all(consistent_silhouette_data$name == seurat_silhouette_data$name))

consistent_silhouette_data$seurat_cluster = seurat_silhouette_data$seurat_cluster

consistent_silhouette_data$closest = ifelse(consistent_silhouette_data$sil_width > 0, consistent_silhouette_data$seurat_cluster, consistent_silhouette_data$neighbor)

consistent_silhouette_data = consistent_silhouette_data |>
  dplyr::select(name, seurat_cluster, closest, sil_width) |>
  dplyr::rename(cluster = seurat_cluster)

write.csv(consistent_silhouette_data, outFile, row.names=FALSE, quote=FALSE)

cat("Saving closest clusters to:", outFile, "\n")  
