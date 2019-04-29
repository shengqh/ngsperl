
rm(list = ls())
setwd("H:/shengquanhu/projects/temp")
source("E:/sqh/programs/perl/ngsperl/lib/Integration/TargetScan7.r")


# miRNA and mRNA files
data <- data.frame(
  smRNA=c(
    "T:/Shared/Labs/Vickers Lab/Tiger/projects/20181205_smallRNA_269_933_2002_human_v4_hg19/host_genome/deseq2_miRNA_TotalReads/result/Cells_24hMetaC_VS_24h5mM_min5_pvalue0.05_DESeq2_sig.csv"),
  mRNA = c(
    "T:/Shared/Labs/Vickers Lab/Tiger/projects/20180822_rnaseq_WZ_371_human/deseq2_proteincoding_genetable/MetaC_vs_5mM_min5_fdr0.05_DESeq2_sig.csv")
, stringsAsFactors = F)
write.csv(data, "Aronoff_mRNA_smRNA_files.csv", row.names = F, quote=F)


targetscan_folder = "T:/Shared/Labs/Vickers Lab/MARS/TargetScan_7_2"

comparison <- c("Cells_24hMetaC_VS_24h5mM")
i<-1
for (i in 1:nrow(data)){
  species = "Hsa"
  prediction = "All"
  project_name = comparison[i]
  miRNA <- fread(data[i,1], data.table = F)
  miRNA <- miRNA[,c("V1", "log2FoldChange")]
  colnames(miRNA) <- c("Feature_miRNA_name", "log2FoldChange_miRNA")
  
  mRNA <- fread(data[i,2], data.table = F)
  mRNA <- mRNA[,c("Feature_gene_name", "log2FoldChange")] 
  colnames(mRNA) <- c("Feature_gene_name", "log2FoldChange_mRNA")
  
  save_targets<-TRUE
  
  targetscan_DE(miRNA, mRNA, targetscan_folder = targetscan_folder, species = species, prediction = prediction, save_targets=save_targets, project_name = project_name)
}
