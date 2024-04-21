
require(biomaRt)
require(stringr)

params_def=read.table(parSampleFile1, stringsAsFactor=F, sep="\t")
params<-split(params_def$V1, params_def$V2)

host=params$host
dataset=params$dataset
symbolKey=params$symbolKey
genesStr<-params$genesStr
addChr<-params$add_chr=="1"

if(!file.exists(genesStr)){
  genesStr = gsub(",", " ", genesStr)
  genesStr = gsub(";", " ", genesStr)
  genes = unlist(strsplit(genesStr, "\\s+"))
}else{
  geneTable = read.table(genesStr, sep="\t", stringsAsFactor=F)
  if (grepl(".bed$", genesStr)){
    genes = unique(geneTable$V4)
    if (all(grepl("(.+)", genes))) {
      genes = gsub(".+?\\(", "", genes)
      genes = gsub("\\).+", "", genes)
      genes = unique(genes)
    }
  }else{
    genes = geneTable$V1
  }
}

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

if(addChr & (!any(grepl("chr", geneLocus$chromosome_name)))){
  geneLocus$chromosome_name = paste0("chr", geneLocus$chromosome_name)
}

bedFile<-paste0(outFile, ".bed")
write.table(geneLocus, file=bedFile, row.names=F, col.names = F, sep="\t", quote=F)

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
