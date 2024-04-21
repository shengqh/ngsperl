
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

cnvtable<-read.table(parFile1, sep="\t", header=T, comment="*", stringsAsFactor=F)
genes<-cnvtable$Gene[!duplicated(cnvtable$Gene)]
genes<-unlist(lapply(genes, function(x) gsub(";.+", "", x)))
genes<-genes[!duplicated(genes)]

ensembl <- useMart("ensembl", host=host, dataset=dataset)

geneLocus<-getBM(attributes=c("chromosome_name", "start_position", "end_position", symbolKey, "strand", "ensembl_gene_id"),
                 filters=symbolKey, values=genes, mart=ensembl, uniqueRows=TRUE,
                 useCache=FALSE)

geneLocus$score<-1000

geneLocus<-geneLocus[,c("chromosome_name", "start_position", "end_position", "score", symbolKey, "strand", "ensembl_gene_id")]
geneLocus<-geneLocus[order(geneLocus$chromosome_name, geneLocus$start_position),]

if(shift != 0){
  geneLocus$start_position[geneLocus$strand == 1] <- geneLocus$start_position[geneLocus$strand == 1] - shift
  geneLocus$end_position[geneLocus$strand == -1] <- geneLocus$end_position[geneLocus$strand == -1] + shift
}

geneLocus$strand[geneLocus$strand == 1]<-"+"
geneLocus$strand[geneLocus$strand == -1]<-"-"

if(addChr & (!any(grepl("chr", geneLocus$chromosome_name)))){
  geneLocus$chromosome_name = paste0("chr", geneLocus$chromosome_name)
}


bedFile<-paste0(outFile, ".bed")
write.table(geneLocus, file=bedFile, row.names=F, col.names = F, sep="\t", quote=F)

missing<-genes[!(genes %in% geneLocus[,symbolKey])]
writeLines(missing, paste0(outFile, ".missing"))

exonLocus<-getBM(attributes=c("exon_chrom_start", "exon_chrom_end", "ensembl_gene_id"),
                 filters="ensembl_gene_id", values=geneLocus$ensembl_gene_id, mart=ensembl, uniqueRows=TRUE,
                 useCache=FALSE)

exonLocus<-exonLocus[,c("ensembl_gene_id", "exon_chrom_start", "exon_chrom_end")]

write.table(exonLocus, file=paste0(bedFile, ".exon"), row.names=F, col.names = F, sep="\t", quote=F)
