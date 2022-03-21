library("MatrixEQTL")

args = commandArgs(trailingOnly=TRUE)

snp_genotype_file=args[1]
snp_location_file=args[2]
gene_expression_file=args[3]
gene_location_file=args[4]
sort_data_by_tcga=args[5]
output_cis_file=args[6]
output_trans_file=args[7]

cat("snp_genotype_file=", snp_genotype_file, "\n")
cat("snp_location_file=", snp_location_file, "\n")
cat("gene_expression_file=", gene_expression_file, "\n")
cat("gene_location_file=", gene_location_file, "\n")
cat("sort_data_by_tcga=", sort_data_by_tcga, "\n")
cat("output_cis_file=", output_cis_file, "\n")
cat("output_trans_file=", output_trans_file, "\n")

pvalue <- 1e-2

checkFileExists<-function(fileName){
  if(!file.exists(fileName)){
    stop(paste0("File not exists ", fileName))
  }
}

checkFileExists(snp_genotype_file)
checkFileExists(snp_location_file)
checkFileExists(gene_expression_file)
checkFileExists(gene_location_file)

useModel = modelLINEAR

if(!is.na(output_trans_file)){
  # perform both cis and trans
  pvOutputThreshold_tra = pvalue
}else{
  # perform local (cis) only
  pvOutputThreshold_tra = 0;
  output_trans_file = ""
}

pvOutputThreshold_cis = pvalue;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

snp_file=snp_genotype_file
exp_file=gene_expression_file
if(sort_data_by_tcga == "1"){
  snpData<-read.table(snp_genotype_file, sep="\t", stringsAsFactor=F, header=T, check.names=F)
  snpNames<-substring(colnames(snpData)[2:ncol(snpData)], 1, 12)
  colnames(snpData)<-c("id", snpNames)

  expData<-read.table(gene_expression_file, sep="\t", stringsAsFactor=F, header=T, check.names=F)
  expNames<-substring(colnames(expData)[2:ncol(expData)], 1, 12)
  colnames(expData)<-c("ID", expNames)
  
  commonNames<-snpNames[snpNames %in% expNames]
  
  snpData<-snpData[,c("id", commonNames)]
  snp_file<-paste0(basename(snp_file), ".tmp.snp")
  write.table(snpData, file=snp_file, quote=F, sep="\t", row.names=F)
  
  expData<-expData[,c("ID", commonNames)]
  exp_file<-paste0(basename(gene_expression_file), ".tmp.geneexpression")
  write.table(expData, file=exp_file, quote=F, sep="\t", row.names=F)
}

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(snp_file);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(exp_file);

## Run the analysis
snpspos = read.table(snp_location_file, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  output_file_name     = output_trans_file,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_cis_file,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

if(sort_data_by_tcga == "1"){
  unlink(snp_file)
  unlink(exp_file)
}
