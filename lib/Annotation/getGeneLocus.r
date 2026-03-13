rm(list=ls()) 
outFile='AutoGVP'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/cindy_chen_projects/20260311_AutoGVP/genes_locus/result')

### Parameter setting end ###

require(biomaRt)
require(stringr)

params_def=read.table(parSampleFile1, stringsAsFactor=F, sep="\t")
params<-split(params_def$V1, params_def$V2)

host=params$host
dataset=params$dataset
symbolKey=params$symbolKey
genesStr=params$genesStr
addChr=params$add_chr=="1"
gene_shift=as.numeric(params$gene_shift)
output_gff=params$output_gff=="1"

if(!file.exists(genesStr)){
  genes = str_trim(unlist(strsplit(genesStr, ",")))
  genesStr = gsub(",", " ", genesStr)
  genesStr = gsub(";", " ", genesStr)
  genes = strsplit(genesStr, "\\s+")
}else{
  geneTable = read.table(genesStr, sep="\t", stringsAsFactor=F)
  if (grepl(".bed$", genesStr)){
    genes = geneTable$V4
    if (all(grepl("(.+)", genes))) {
      genes = gsub(".+?\\(", "", genes)
      genes = gsub("\\).+", "", genes)
    }
  }else{
    genes = geneTable$V1
  }
}
genes = genes |> unlist() |> str_trim() |> unique() |> sort()

cat("genes:", paste0(genes, collapse=", "), "\n")

ensembl <- useMart("ensembl", host=host, dataset=dataset)

geneLocus<-getBM(attributes=c("chromosome_name", "start_position", "end_position", symbolKey, "strand", "ensembl_gene_id"),
                 filters=symbolKey, 
                 values=genes, 
                 mart=ensembl, 
                 uniqueRows=TRUE,
                 useCache=FALSE)

geneLocus<-geneLocus[nchar(geneLocus$chromosome_name) < 6,]

geneLocus$score<-1000

geneLocus<-geneLocus[,c("chromosome_name", "start_position", "end_position", "score", symbolKey, "strand", "ensembl_gene_id")]
geneLocus<-geneLocus[order(geneLocus$chromosome_name, geneLocus$start_position),]

geneLocus$strand[geneLocus$strand == 1]<-"+"
geneLocus$strand[geneLocus$strand == -1]<-"-"

if(gene_shift != 0){
  cat("shift gene ", gene_shift, " bases\n")
  geneLocus$start_position = geneLocus$start_position - gene_shift
  geneLocus$start_position[geneLocus$start_position < 1] = 1
  geneLocus$end_position = geneLocus$end_position + gene_shift
}

if(addChr & (!any(grepl("chr", geneLocus$chromosome_name)))){
  geneLocus$chromosome_name = paste0("chr", geneLocus$chromosome_name)
}

geneLocus$chromosome_name = gsub("^chrMT$", "chrM", geneLocus$chromosome_name)
geneLocus$chromosome_name = gsub("^MT$", "M", geneLocus$chromosome_name)

bedFile<-paste0(outFile, ".bed")
geneLocusBed = geneLocus |> dplyr::mutate(start_position = start_position-1)
write.table(geneLocusBed, file=bedFile, row.names=F, col.names = F, sep="\t", quote=F)

missing<-genes[!(genes %in% geneLocus[,symbolKey])]
writeLines(missing, paste0(outFile, ".missing"))

exonLocus<-getBM(attributes=c("exon_chrom_start", "exon_chrom_end", "ensembl_gene_id"),
                 filters="ensembl_gene_id", 
                 values=geneLocus$ensembl_gene_id, 
                 mart=ensembl, 
                 uniqueRows=TRUE,
                 useCache=FALSE)

exonLocus<-exonLocus[,c("ensembl_gene_id", "exon_chrom_start", "exon_chrom_end")]

write.table(exonLocus, file=paste0(bedFile, ".exon"), row.names=F, col.names = F, sep="\t", quote=F)

if(output_gff){
  geneLocus$unknown1="GENE"
  geneLocus$unknown2='.'
  geneLocus$unknown3='.'
  geneLocus=geneLocus[,c(1,4,6,2,3,7,5,8)]
  write.table(geneLocus, file=paste0(outFile, ".gff"), row.names = F,col.names = FALSE, quote = F, sep="\t")
}
