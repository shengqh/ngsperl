rm(list=ls()) 
outFile='5057_AD_combined'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/workspace/shengq2/kasey_vickers_projects/2020_projects/20200805_5057_AD_rnaseq_hsammu_combined_byMars/20240401_bamplot/bamplot_gene_gff/result')

### Parameter setting end ###

source("countTableVisFunctions.R")
options(bitmapType='cairo')

library(biomaRt)
library(gtools)
library(stringr)


myoptions<-read_file_map(parSampleFile1, do_unlist=FALSE)

outFile = myoptions$task_name

biomart_host = myoptions$biomart_host
if(!grepl("^http", biomart_host)){
  biomart_host = paste0("https://", biomart_host)
}
biomart_dataset = myoptions$biomart_dataset
biomart_symbolKey = myoptions$biomart_symbolKey
biomart_add_chr = is_one(myoptions$biomart_add_chr)
biomart_add_prefix = myoptions$biomart_add_prefix
prefix=paste0(biomart_add_prefix, ifelse(biomart_add_chr, "chr", ""))

output_gff = is_one(myoptions$output_gff)

gene_names = myoptions$gene_names
genes = str_trim(unlist(strsplit(gene_names, ","))) |> sort()
gene_shift = as.numeric(myoptions$gene_shift)

cat("gene_names:", gene_names, "\n")

cat("connect to", biomart_host, "for", biomart_dataset, "\n")
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host=biomart_host, dataset=biomart_dataset)
allgenepos <- lapply(genes, function(x){
  getBM(attributes = c("chromosome_name", "start_position", "end_position", biomart_symbolKey, "strand"),
        filters=biomart_symbolKey,
        values = x,
        mart=mart)})

gene_locus <- do.call("rbind", allgenepos)
rownames(gene_locus)<-gene_locus$ID
colnames(gene_locus)<-c("chr",  "start",  "end", "gene", "strand")
gene_locus$strand[gene_locus$strand==1] = '+'
gene_locus$strand[gene_locus$strand==-1] = '-'
  
gene_locus$chr = paste0(prefix, gene_locus$chr)

write.table(gene_locus, file=paste0(myoptions$task_name, ".gene_locus_raw.txt"), row.names = F,col.names = TRUE, quote = F, sep="\t")

if(gene_shift != 0){
  cat("shift gene ", gene_shift, " bases\n")
  is_forward_strand= gene_locus$strand == '+'
  gene_locus$start = gene_locus$start - gene_shift
  gene_locus$start[gene_locus$start < 1] = 1
  gene_locus$end = gene_locus$end + gene_shift
}

ll = unlist(lapply(gene_locus$chr, nchar))
gene_locus = gene_locus[ll < 6,]

missed_file=paste0(myoptions$task_name, ".missed_genes.txt")
if(file.exists(missed_file)){
  unlink(missed_file)
}

missed_genes = genes[!genes %in% gene_locus$gene]
if(length(missed_genes) > 0){
  cat("missed genes:", missed_genes, "\n")
  writeLines(missed_genes, missed_file)
}

if(output_gff){
  gene_locus$unknown1="GENE"
  gene_locus$unknown2='.'
  gene_locus$unknown3='.'
  gene_locus=gene_locus[,c(1,4,6,2,3,7,5,8)]
  outFile = paste0(myoptions$task_name, ".gff")
  col.names = FALSE
}else{
  outFile = paste0(myoptions$task_name, ".txt")
  col.names = TRUE
}
cat("write to", outFile, "\n")
write.table(gene_locus, file=outFile, row.names = F,col.names = col.names, quote = F, sep="\t")
