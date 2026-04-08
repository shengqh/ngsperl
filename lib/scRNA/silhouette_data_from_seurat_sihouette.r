args <- commandArgs(trailingOnly = TRUE)

# args=c("/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_fastmnn_dr0.5_3_choose_silhouette/result/Aorta_Progeria.seurat.silhouette.csv", 
# "/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_fastmnn_dr0.5_3_choose_silhouette/result/Aorta_Progeria.seurat.coords.csv",
# "Aorta_Progeria")

if (length(args) != 3) {
  stop("Usage: Rscript silhouette_data_from_sihouette.r <silhouette_file> <coords_file> <output_cluster_file>")
}

silhouette_file <- args[1]
coords_file <- args[2]
output_cluster_file <- args[3]

cat("silhouette_file:", silhouette_file, "\n")
cat("coords_file:", coords_file, "\n")
cat("output_cluster_file:", output_cluster_file, "\n")

library(data.table)

cat("Reading silhouette file:", silhouette_file, "\n")
silhouette_data <- fread(silhouette_file)

cat("Reading coords file:", coords_file, "\n")
coords_data <- data.frame(fread(coords_file), row.names=1)

stopifnot(all(rownames(coords_data) %in% silhouette_data$name))

silhouette_data$new_cluster = ifelse(silhouette_data$sil_width > 0, silhouette_data$cluster, -1)

writeLines(as.character(silhouette_data$new_cluster), output_cluster_file)

cat("Saving consistent clusters to:", output_cluster_file, "\n")  
