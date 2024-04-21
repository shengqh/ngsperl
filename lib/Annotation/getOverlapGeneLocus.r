
require(biomaRt)
require(stringr)

if(!exists("host")){
  host="grch37.ensembl.org"
}

if(!exists("dataset")){
  dataset = "hsapiens_gene_ensembl"
}

if(!exists("symbolKey")){
  symbolKey = "hgnc_symbol"
}

if(!exists('shift')){
  shift = 0
}

if(!exists('outFile')){
  outFile='test'
}

if(!exists('addChr')){
  addChr = 0
}

ensembl <- useMart("ensembl", host=host, dataset=dataset)

coords = read.table(parFile1, sep="\t", header=T, stringsAsFactors = F)
filterlist=coords$Locus
filterlist<-gsub("-",":",filterlist)
filterlist<-gsub("^chr","",filterlist)

geneLocus=getBM(attributes = c("chromosome_name", "start_position", "end_position", symbolKey, "strand", "ensembl_gene_id"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), 
              mart = ensembl,
              useCache=FALSE)

geneLocus$score<-1000

geneLocus<-geneLocus[,c("chromosome_name", "start_position", "end_position", "score", symbolKey, "strand", "ensembl_gene_id")]
geneLocus<-geneLocus[order(geneLocus$chromosome_name, geneLocus$start_position),]
geneLocus$strand[geneLocus$strand == 1]<-"+"
geneLocus$strand[geneLocus$strand == -1]<-"-"

if(addChr & (!any(grepl("chr", geneLocus$chromosome_name)))){
  geneLocus$chromosome_name = paste0("chr", geneLocus$chromosome_name)
}

bedFile<-paste0(outFile, ".bed")
write.table(geneLocus, file=bedFile, row.names=F, col.names = F, sep="\t", quote=F)

exonLocus<-getBM(attributes=c("exon_chrom_start", "exon_chrom_end", "ensembl_gene_id"),
                 filters="ensembl_gene_id", 
                 values=geneLocus$ensembl_gene_id, 
                 mart=ensembl, 
                 uniqueRows=TRUE,
                 useCache=FALSE)

exonLocus<-exonLocus[,c("ensembl_gene_id", "exon_chrom_start", "exon_chrom_end")]
write.table(exonLocus, file=paste0(bedFile, ".exon"), row.names=F, col.names = F, sep="\t", quote=F)
