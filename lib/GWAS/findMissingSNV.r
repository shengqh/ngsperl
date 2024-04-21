require(biomaRt)
require(stringr)

args = commandArgs(trailingOnly=TRUE)
riskDbsnpFile = args[1]
plinkBimFile = args[2]

outputfile = args[2]


#setwd("H:/shengquanhu/projects/20190320_2118-JB_GSProject_GRS")
setwd("/data/h_vangard_1/macrae_linton_data/2118")

ensembl <- useMart("ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset = "hsapiens_snp")

riskDbsnp<-read.table(riskDbsnpFile, sep="\t")
riskDbsnp<-riskDbsnp[riskDbsnp$V1 != "",]
rownames(riskDbsnp)<-riskDbsnp$V1

mappedDbsnp<-getBM(attributes=c(
  "refsnp_id", "chr_name", "chrom_start", "chrom_end", "chrom_strand",
  "allele", "allele_1", "minor_allele"),
  filters="snp_filter", 
  values=riskDbsnp$V1,
  mart=ensembl, 
  uniqueRows=TRUE,
  useCache=FALSE)
mappedDbsnp<-mappedDbsnp[str_length(mappedDbsnp$chr_name) < 3,]
rownames(mappedDbsnp)<-mappedDbsnp$refsnp_id
mappedDbsnp<-mappedDbsnp[rownames(riskDbsnp),]
mappedDbsnp$locus<-paste0(mappedDbsnp$chr_name,":",mappedDbsnp$chrom_start)

merged<-merge(x = riskDbsnp, y = mappedDbsnp, by.x = "V1", by.y="refsnp_id", all.x = TRUE)

illumina<-read.table(plinkBimFile, sep="\t")
illumina$locus<-paste0(illumina$V1,":",illumina$V4)

illuminaMatch<-illumina[illumina$locus %in% mappedDbsnp$locus,]

multiIllumina<-table(illuminaMatch$locus)[table(illuminaMatch$locus) > 1]
illuminaMatchUnique<-illuminaMatch[!(illuminaMatch$locus %in% names(multiIllumina)),]
illuminaMatchMulti<-illuminaMatch[(illuminaMatch$locus %in% names(multiIllumina)),]
illuminaMatchMultiUniq<-illuminaMatchMulti[grepl("^rs", illuminaMatchMulti$V2),]

finalIlluminaMatch<-rbind(illuminaMatchUnique, illuminaMatchMultiUniq)


mergedAll<-merge(y = finalIlluminaMatch, x= merged, by = "locus", all.x = TRUE)
mergedMatch<-mergedAll[!is.na(mergedAll$V1.y),]
finalTable<-mergedMatch[,c("V2.y", "V2.x", "V3.x")]

finalTable<-finalTable[order(finalTable$V2.y),]
write.table(finalTable, file="LDL_GLGC81SNPs_list_effectsize.update.raw", sep="\t", row.names=F, col.names=F, quote = F)


imputationTable<-mergedAll[is.na(mergedAll$V1.y), ]
imputationTable<-imputationTable[,c("chr_name", "chrom_start", "chrom_end", "V1.x")]
imputationTable$chr_name<-as.numeric(imputationTable$chr_name)
imputationTable$chrom_start<-as.numeric(imputationTable$chrom_start)
imputationTable<-imputationTable[order(imputationTable$chr_name, imputationTable$chrom_start),]

write.table(imputationTable, file="LDL_GLGC81SNPs_list_effectsize.missing.bed", sep="\t", row.names=F, col.names=F, quote = F)
