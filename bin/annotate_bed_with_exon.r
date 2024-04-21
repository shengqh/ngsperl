library(biomaRt)
library(GenomicRanges)
library(rtracklayer)

#setwd("H:/shengquanhu/projects/database/exomeseq/Twist")
setwd("/scratch/cqs/references/exomeseq/Twist")

beddata<-read.table("Twist_Exome_Target_hg38.slop50.bed", sep="\t", stringsAsFactors=F)

host="https://www.ensembl.org";dataset="hsapiens_gene_ensembl";symbolKey="hgnc_symbol";
ensembl <- useMart("ensembl", host=host, dataset=dataset)
exons <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand"), mart = ensembl)
exons$chromosome_name=paste0("chr", exons$chromosome_name)

genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = ensembl)
dg<-genes[duplicated(genes$ensembl_gene_id),]
nodg<-genes[!duplicated(genes$ensembl_gene_id),]
nodg_names<-split(nodg$hgnc_symbol, nodg$ensembl_gene_id)

hg38.exons<-GRanges(exons$chromosome_name, ranges=IRanges(start=exons$exon_chrom_start, end=exons$exon_chrom_end, names = exons$ensembl_exon_id))
hg39.exome<-GRanges(beddata$V1, ranges=IRanges(start=beddata$V2, end=beddata$V3))

ov <- findOverlaps(hg39.exome, hg38.exons, type = "any", select="first")
ov.na<-which(is.na(ov))
ov.not.na<-which(!is.na(ov))

exome_exon<-beddata
exome_exon$name<-paste0(exome_exon$V1, "_", exome_exon$V2, "_", exome_exon$V3)
exome_exon$name[ov.not.na] = paste0(exome_exon$name[ov.not.na], "(", nodg_names[exons$ensembl_gene_id[ov[ov.not.na]]], ")")
exome_exon$score<-'0'
exome_exon$strand<-'*'
exome_exon$strand[ov.not.na] = exons$strand[ov[ov.not.na]]
exome_exon$strand[exome_exon$strand=="1"] = '+'
exome_exon$strand[exome_exon$strand=="-1"] = '-'

write.table(exome_exon, file="Twist_Exome_Target_hg38.slop50.name.bed", sep="\t", row.names=F, col.names=F, quote=F)

exonsaf<-data.frame("GeneID"=exome_exon$name,
                    "Chr"=exome_exon$V1,
                    "Start"=exome_exon$V2 + 1,
                    "End"=exome_exon$V3,
                    "Strand"=exome_exon$strand)
write.table(exonsaf, file="Twist_Exome_Target_hg38.slop50.name.saf", sep="\t", row.names=F, col.names=T, quote=F)

