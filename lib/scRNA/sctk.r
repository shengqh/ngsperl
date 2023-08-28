rm(list=ls()) 
sample_name='Adipose_9240'
outFile='Adipose_9240'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230410_combined_scRNA_hg38_test_vst2/decontX_sctk/result/Adipose_9240')

### Parameter setting end ###

source("scRNA_func.r")
library("singleCellTK")
library("SingleCellExperiment")
library("Seurat")

sample_map = read_file_map("fileList1.txt")

myoptions = read_file_map("fileList2.txt", do_unlist=FALSE)
nFeature_cutoff_min=max(as.numeric(myoptions$nFeature_cutoff_min) - 50, 100)
nCount_cutoff=max(as.numeric(myoptions$nCount_cutoff) - 50, 200)

overwrite=FALSE

sce_list <- list()
sample_name=names(sample_map)[1]
for (sample_name in names(sample_map)){
  cat("sample_name\n")
  countfile = sample_map[sample_name]

  sctk_rds = paste0(sample_name, ".sctk.rds")

  redo = !file.exists(sctk_rds) | overwrite
  if(!redo){
    cat("  read", sctk_rds, "\n")
    sce = readRDS(sctk_rds)
  }else{
    cat("  read_scrna_data\n")
    lst = read_scrna_data(countfile)
    counts<-lst$counts
    rawobj = CreateSeuratObject(counts = counts, project = sample_name)
    rawobj<-subset(rawobj, subset = nFeature_RNA >= nFeature_cutoff_min & nCount_RNA >= nCount_cutoff)

    sce <- as.SingleCellExperiment(rawobj)
    rm(rawobj)

    cat("  runCellQC\n")
    #sce = runCellQC(sce, algorithms= c("QCMetrics", "scrublet", "doubletFinder", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid", "decontX"))
    #There is error to run scrublet, so ignore it.
    sce = runCellQC(sce, algorithms= c("QCMetrics", "doubletFinder", "scDblFinder", "cxds", "bcds", "cxds_bcds_hybrid", "decontX"))

    saveRDS(sce, sctk_rds)
  }

  saveRDS(colData(sce), paste0(sample_name, ".meta.rds"))

  #cat("  reportCellQC\n")
  #there is a lot of error to run reportCellQC, so ignore it.
  #reportCellQC(inSCE = sce, output_file = paste0(sample_name, "_reportCellQC"))

  #cat("  exportSCE\n")
  #exportSCE(inSCE = sce, samplename = sample_name, type = "Cells", format = c("Seurat"))
}
