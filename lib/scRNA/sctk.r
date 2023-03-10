library("singleCellTK")
library("SingleCellExperiment")
source("scRNA_func.r")

sample_dirs <- read.table("fileList1.txt", header=F)
sample_map = unlist(split(sample_dirs$V1, sample_dirs$V2))

sce_list <- list()
sample_name=names(sample_map)[1]
for (sample_name in names(sample_map)){
  cat("sample_name\n")
  countfile = sample_map[sample_name]

  cat("  read_scrna_data\n")
  lst = read_scrna_data(countfile)
  counts<-lst$counts
  sce <- SingleCellExperiment(assays=list(counts=counts))

  cat("  runCellQC\n")
  sce = runCellQC(sce, algorithms= c("QCMetrics", "scrublet", "doubletFinder", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid", "decontX"))

  saveRDS(sce, paste0(sample_name, ".sctk.rds"))

  cat("  reportCellQC\n")
  reportCellQC(inSCE = sce, output_file = paste0(sample_name, "_reportCellQC"))

  cat("  exportSCE\n")
  exportSCE(inSCE = sce, samplename = sample_name, type = "Cells", format = c("Seurat"))
}
