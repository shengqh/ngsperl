library(rsnps)
library(dplyr)

setwd("C:/projects/jonathan_brown/20200905_wgs_5162_hg38/")

snps=c("rs738409", "rs58542926", "rs641738", "rs72613567", "rs738408", "rs3747207", "rs5748926", "rs698718", "rs1260326", "rs61756425")

gts<-lapply(snps, function(x){
  ncbi_snp_query(x)
})

gtf=bind_rows(gts)

gtf2=gtf[,c(2,3,3,5)]
gtf2$chromosome=as.numeric(gtf2$chromosome)
gtf2$bp=as.numeric(gtf2$bp)
gtf2$bp=gtf2$bp-1
gtf2=gtf2[order(gtf2$chromosome, gtf2$bp),]

gtf2$chromosome = paste0("chr", gtf2$chromosome)
write.table(gtf2, "snps.bed", sep="\t", row.names=F, col.names=F, quote=F)

 