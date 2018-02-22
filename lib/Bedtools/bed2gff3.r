library(rtracklayer)

args <- commandArgs(TRUE)

print(args)

bedfile=args[1]
gff3file=args[2]

## import the bed file
bed.ranges <- import.bed(bedfile)

## export as a gff3 file
export.gff3(bed.ranges, gff3file)

gf<-read.table(gff3file, sep="\t", header=F)
gf$V2<-paste0("Peak", c(1:nrow(gf)))

write.table(gf, file=gff3file, sep="\t", col.names=F, row.names=F, quote=F)
