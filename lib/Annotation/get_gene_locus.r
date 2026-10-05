rm(list=ls()) 
outFile='VectorSeq'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20260713_MVendo_RNAseq_hg38_vector/genes_locus/result')

### Parameter setting end ###

library(AnnotationHub)
library(ensembldb)
library(stringr)

if(!dir.exists("./AnnotationHub_cache")){
  dir.create("./AnnotationHub_cache")
}

setAnnotationHubOption("CACHE", "./AnnotationHub_cache")

params_def=read.table(parSampleFile1, stringsAsFactor=F, sep="\t")
params<-split(params_def$V1, params_def$V2)
genesStr=params$genesStr
frank_bases=as.numeric(params$gene_shift)
output_gff=params$output_gff=="1"

gene_names <- trimws(strsplit(genesStr, ",")[[1]])

cat("gene_names: ", gene_names, "\n")

db_str="EnsDb:Homo sapiens:113"
db <- trimws(strsplit(db_str, ":")[[1]])

cat("db: ", db, "\n")

addChr=params$add_chr=="1"

ah <- AnnotationHub()

edb = query(ah, db)
edb <- edb[[1]]

geneLocus = genes(
    edb,
    filter = AnnotationFilterList(GeneNameFilter(gene_names),
                                  GeneBiotypeFilter("protein_coding")),
    return.type = "DataFrame"
)

#          gene_id   gene_name   gene_biotype gene_seq_start gene_seq_end
#       <character> <character>    <character>      <integer>    <integer>
# 1 ENSG00000149311         ATM protein_coding      108222804    108369102
#      seq_name seq_strand seq_coord_system            description
#   <character>  <integer>      <character>            <character>
# 1          11          1       chromosome ATM serine/threonine..
#      gene_id_version canonical_transcript      symbol entrezid
#          <character>          <character> <character>   <list>
# 1 ENSG00000149311.22      ENST00000675843         ATM      472

geneLocus<-geneLocus[nchar(geneLocus$seq_name) < 6,]

geneLocus$score<-1000

geneLocus<-geneLocus[,c("seq_name", "gene_seq_start", "gene_seq_end", "score", "symbol", "seq_strand", "gene_id")]
geneLocus<-geneLocus[order(geneLocus$seq_name, geneLocus$gene_seq_start),]

geneLocus$seq_strand[geneLocus$seq_strand == 1]<-"+"
geneLocus$seq_strand[geneLocus$seq_strand == -1]<-"-"

if(addChr & (!any(grepl("chr", geneLocus$seq_name)))){
  geneLocus$seq_name = paste0("chr", geneLocus$seq_name)
}

geneLocus$seq_name=gsub("chrMT", "chrM", geneLocus$seq_name)

if(frank_bases > 0){
  geneLocus$gene_seq_start = geneLocus$gene_seq_start - frank_bases
  geneLocus$gene_seq_end = geneLocus$gene_seq_end + frank_bases
}

bedFile<-paste0(outFile, ".bed")
write.table(geneLocus, file=bedFile, row.names=F, col.names = F, sep="\t", quote=F)

interval<-paste0(geneLocus$seq_name[1], ":", geneLocus$gene_seq_start[1], "-", geneLocus$gene_seq_end[1])
writeLines(interval, con=paste0(outFile, ".interval"))

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
