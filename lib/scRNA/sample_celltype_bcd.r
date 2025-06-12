rm(list=ls()) 
sample_name='CD_Met_01'
outFile='CD_Met_01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/nd_seurat_sct2_merge_dr0.2_3_choose/result/P10940.meta.rds'
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/20250606_homer/T01_sample_celltype_bcd/result/CD_Met_01')

### Parameter setting end ###

source("scRNA_func.r")
source("reportFunctions.R")
myoptions = read_file_map(parSampleFile2, do_unlist=FALSE)

group_by_column = myoptions$group_by_column
cat("group_by_column =", group_by_column, "\n")

sample_column = myoptions$sample_column
cat("sample_column =", sample_column, "\n")

original_cell_column = myoptions$original_cell_column
cat("original_cell_column =", original_cell_column, "\n")

meta=readRDS(parFile1)

if (sample_column %in% colnames(meta) == FALSE) {
  stop(paste("sample_column", sample_column, "not found in metadata!"))
}
samples = unique(meta[[sample_column]])

if(!sample_name %in% samples) {
  stop(paste("sample_name", sample_name, "not found in metadata!"))
}

if(group_by_column %in% colnames(meta) == FALSE) {
  stop(paste("group_by_column", group_by_column, "not found in metadata!"))
}

cat("Processing sample:", sample_name, "\n")
meta_sample = meta[meta[[sample_column]] == sample_name, ]
if (original_cell_column == "rownames"){
  bcd=data.frame(cell=rownames(meta_sample), celltype=meta_sample[[celltype_column]])
}else{
  if (original_cell_column %in% colnames(meta_sample) == FALSE) {
    stop(paste("original_cell_column", original_cell_column, "not found in metadata!"))
  }
  bcd=data.frame(cell=meta_sample[[original_cell_column]], celltype=paste0(sample_name, ".", celltype_to_filename(meta_sample[[celltype_column]])))
}

sample_file = paste0(sample_name, ".bcd.txt")
write.table(bcd, file=sample_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

