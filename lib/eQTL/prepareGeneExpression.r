library(preprocessCore)
library(GenABEL)
library(gtools)

args = commandArgs(trailingOnly=TRUE)

DEBUG<-F
if(DEBUG){
  #fpkmfile<-"/workspace/shengq1/tcga/brca/normal_table/result/brca_rnaseq_normal_fpkm.tsv"
  #commfile<-"/workspace/shengq1/tcga/brca/eqtl/data/normal_brca_common.fam"
  #posfile<-"/workspace/shengq1/tcga/hg19.pos"
  #outputfile<-"/workspace/shengq1/tcga/brca/eqtl/fastqtl/normal_brca_gene.bed"
  fpkmfile<-"/scratch/cqs/shengq1/guoyan/prepareRnaseq/result/BRCA/BRCA.rnaseq2.fpkm.normal.tsv"
  commfile<-"/scratch/cqs/shengq1/guoyan/dataPreparationNormal/result/BRCA/normal_BRCA_common.fam"
  posfile<-"/scratch/cqs/shengq1/guoyan/prepareRnaseq/result/BRCA/BRCA.rnaseq2.gene.pos"
  outputfile<-"/scratch/cqs/shengq1/guoyan/dataPreparationNormal/result/BRCA/normal_BRCA_gene.bed"
}else{
  fpkmfile=args[1]
  commfile=args[2]
  posfile<-args[3]
  outputfile=args[4]
}

cat("fpkmfile=", fpkmfile, "\n")
cat("commfile=", commfile, "\n")
cat("posfile=", posfile, "\n")
cat("outputfile=", outputfile, "\n")

#http://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods

data<-read.table(fpkmfile, stringsAsFactors = F, sep="\t", header=T, row.names=1, check.names = F)
samples<-read.table(commfile,sep=" ",stringsAsFactors = F)

fpkmdata<-data[,samples$V1]

#GTEx: Genes were selected based on expression thresholds of >0.1 RPKM in at least 10 individuals and ≥6 reads in at least 10 individuals
#Since we don't have raw count, we use only first criteria
fpkmdata.num01<-apply(fpkmdata, 1, function(x){
  sum(x > 0.1)
})
fpkmdata<-fpkmdata[fpkmdata.num01 >= 10,]

if(nrow(fpkmdata) == 0){
  stop("No gene passed the criteria!");
}

#Expression values were quantile normalized to the average empirical distribution observed across samples.
fpkmdata.nq=normalize.quantiles(as.matrix(fpkmdata))

#For each gene, expression values were inverse quantile normalized to a standard normal distribution across samples.
fpkmdata.int=apply(fpkmdata.nq, 1, rntransform)

fpkmdata.intfm<-data.frame(t(fpkmdata.int))
colnames(fpkmdata.intfm)<-colnames(fpkmdata)
rownames(fpkmdata.intfm)<-rownames(fpkmdata)

allgenepos = read.table(posfile, header=T, sep="\t")
if ("geneid" %in% colnames(allgenepos)){
  colnames(allgenepos)<-c("ID", "chr", "start", "end")
  allgenepos<-allgenepos[,c("chr", "start", "end", "ID")]
}
rownames(allgenepos)<-allgenepos$ID

commonGenes<-as.character(allgenepos$ID[allgenepos$ID %in% rownames(fpkmdata.intfm)])
commonPos<-allgenepos[commonGenes,]
commonData<-fpkmdata.intfm[commonGenes,]
finaldata<-cbind(commonPos, commonData)

finaldata<-finaldata[order(finaldata$chr, finaldata$start ),]
finaldata<-finaldata[mixedorder(finaldata$chr),]
colnames(finaldata)[1]="#Chr"

write.table(finaldata, file=outputfile, row.names = F,col.names = T,  quote = F, sep="\t")
