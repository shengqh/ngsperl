args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript silhouette_pairwise_data_from_seurat.r <seurat_rds> <cluster_column> <reduction_name> <output_prefix>")
}


obj_file <- args[1]
cluster <- args[2]
reduction <- args[3]
outFile <- args[4]

if(0){
  obj_file='/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/cellbender_nd_seurat_fastmnn_dr0.5_3_choose/result/Aorta_Progeria.final.rds'
  cluster='seurat_clusters'
  reduction='umap'
  outFile='Aorta_Progeria.umap.seurat'
}

cat("obj_file:", obj_file, "\n")
cat("cluster:", cluster, "\n")
cat("reduction:", reduction, "\n")
cat("outFile:", outFile, "\n")

library(Seurat)

cat("Reading object file:", obj_file, "\n")
obj <- readRDS(obj_file)

clusters=FetchData(obj, vars = cluster)[, 1]

coords <- Embeddings(obj, reduction = reduction)

writeLines(as.character(clusters), paste0(outFile, ".clusters.txt"))
write.csv(coords, paste0(outFile, ".coords.csv"), row.names=TRUE, quote=FALSE)

cat("Saving clusters and coords to:", paste0(outFile, ".clusters.txt"), "and", paste0(outFile, ".coords.csv"), "\n")  
