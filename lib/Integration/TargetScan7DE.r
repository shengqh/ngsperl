
data<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactor=F)

targetscan_folder = parFile1

comparisons <- unique(data$V3)

comparison<-comparisons[1]
for ( comparison in comparisons ) {
  cdata<-data[data$V3 == comparison,]
  miRNA_file = cdata[cdata$V2=="miRNA", "V1"]
  mRNA_file = cdata[cdata$V2=="mRNA", "V1"]
  species = cdata[cdata$V2=="species", "V1"]

  prediction = "All"
  project_name = comparison
  miRNA <- fread(miRNA_file, data.table = F)
  miRNA <- miRNA[,c("V1", "log2FoldChange")]
  colnames(miRNA) <- c("Feature_miRNA_name", "log2FoldChange_miRNA")

  mRNA <- fread(mRNA_file, data.table = F)
  mRNA <- mRNA[,c("Feature_gene_name", "log2FoldChange")]
  colnames(mRNA) <- c("Feature_gene_name", "log2FoldChange_mRNA")

  save_targets<-TRUE

  targetscan_DE(miRNA, mRNA, targetscan_folder = targetscan_folder, species = species, prediction = prediction, save_targets=save_targets, project_name = project_name)
}

