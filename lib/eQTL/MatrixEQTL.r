library("MatrixEQTL")

args = commandArgs(trailingOnly=TRUE)

snp_genotype_file=args[1]
snp_location_file=args[2]
gene_expression_file=args[3]
gene_location_file=args[4]
output_cis_file=args[5]
output_trans_file=args[6]

cat("snp_genotype_file=", snp_genotype_file, "\n")
cat("snp_location_file=", snp_location_file, "\n")
cat("gene_expression_file=", gene_expression_file, "\n")
cat("gene_location_file=", gene_location_file, "\n")
cat("output_cis_file=", output_cis_file, "\n")
cat("output_trans_file=", output_trans_file, "\n")

useModel = modelLINEAR

# perform local (cis) only
pvOutputThreshold_tra = 0;
pvOutputThreshold_cis = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(snp_genotype_file);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(gene_expression_file);

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
