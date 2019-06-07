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

if(!exists('genesStr')){
  genesStr<-"LDLR APOB PCSK9 LDLRAP1 STAP1 LIPA ABCG5 ABCGB APOE LPA PNPLA5 CH25H INSIG2 SIRT1"
}

if(!exists('outFile')){
  outFile='test'
}

if(!file.exists(genesStr)){
  genes=unlist(strsplit(genesStr, "\\s+"))
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

geneLocus<-getBM(attributes=c("chromosome_name", "start_position", "end_position", symbolKey, "strand"),
                 filters=symbolKey, values=genes, mart=ensembl, uniqueRows=TRUE)

geneLocus$score<-1000
geneLocus<-geneLocus[,c("chromosome_name", "start_position", "end_position", "score", symbolKey, "strand")]
geneLocus<-geneLocus[order(geneLocus$chromosome_name, geneLocus$start_position),]
geneLocus$strand[geneLocus$strand == 1]<-"+"
geneLocus$strand[geneLocus$strand == -1]<-"-"

write.table(geneLocus, file=paste0(outFile, ".bed"), row.names=F, col.names = F, sep="\t", quote=F)

missing<-genes[!(genes %in% geneLocus[,symbolKey])]
writeLines(missing, paste0(outFile, ".missing"))
