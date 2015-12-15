library(DNAcopy)
cn <- read.table(inputfile, header=T, stringsAsFactors = F)
CNA.object <-CNA( genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_start, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=0, min.width=2)
segs2 = segs$output
write.table(segs2[,2:6], file=outputfile, row.names=F, col.names=F, quote=F, sep="\t")