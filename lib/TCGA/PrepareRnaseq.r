library(TCGAbiolinks)
library(SummarizedExperiment)

DEBUG=0

writeFile<-function(dat, file){
  newDat<-cbind(data.frame(ID=rownames(dat)), dat)
  rownames(newDat)<-rownames(dat)
  write.table(newDat, file=file, sep="\t", row.names=F, quote=F)
}

args = commandArgs(trailingOnly=TRUE)
cancerName = args[1]
outputFolder = args[2]

if(DEBUG){
  cancerName = "SKCM"
  outputFolder = "/scratch/cqs/shengq1/guoyan/prepareRnaseq/result/SKCM"
}

setwd(outputFolder)

project = paste0("TCGA-", cancerName)
countFile<-paste0(cancerName, ".rnaseq2.count.rdata")

if(!file.exists(countFile)){
  query <- GDCquery(project = project, 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    experimental.strategy = "RNA-Seq",
                    platform = "Illumina HiSeq",
                    file.type  = "results",
                    legacy = TRUE)
  
  source <- ifelse(query$legacy, "legacy", "harmonized")
  files <- file.path(outputFolder, "GDCdata", query$project, source, gsub(" ", "_", query$results[[1]]$data_category), gsub(" ", "_", query$results[[1]]$data_type), 
                     gsub(" ", "_", query$results[[1]]$file_id), gsub(" ", "_", query$results[[1]]$file_name))
  tcgaNames<-substring(query$results[[1]]$cases,1,15)
  
  tcgaFiles<-data.frame(Name=tcgaNames, File=files)
  write.table(tcgaFiles, paste0(cancerName, ".rnaseq2.filelist"), sep="\t", quote=F, row.names=F)
  
  GDCdownload(query, method = "client", chunks.per.download = 10)
  rnaseq <- GDCprepare(query)
  
  save(rnaseq, file=countFile)
}else{
  load(countFile)
}

pos<-rowRanges(rnaseq)
posframe<-data.frame(pos)
posframe<-posframe[,c("ensembl_gene_id", "seqnames", "start", "end")]
colnames(posframe)<-c("geneid", "chr", "s1", "s2")

exp<-assay(rnaseq)
rownames(exp)<-pos$ensembl_gene_id

isdup<-duplicated(posframe$geneid)
posframe<-posframe[!isdup,]
exp<-exp[!isdup,]

#sort pos by chromosome and start
posframe$chr<-substring(posframe$chr, 4)
posframe$chrX<-grepl('[XYZ]', posframe$chr)
posframe<-posframe[with(posframe, order(chrX, chr, s1, s2)), c("geneid", "chr", "s1", "s2")]

#save position file
write.table(posframe, file=paste0(cancerName, ".rnaseq2.gene.pos"), sep="\t", row.names=F, quote=F)

exp<-exp[as.character(posframe$geneid),]
colnames(exp)<-substring(colnames(exp), 1, 15)

writeFile(exp, file=paste0(cancerName, ".rnaseq2.count.tsv"))

normal<-exp[,grepl("11$", colnames(exp))]
if(is.matrix(normal) && ncol(normal) > 0){
  writeFile(normal, paste0(cancerName, ".rnaseq2.count.normal.tsv"))
}

tumor<-exp[,grepl("01$", colnames(exp))]
if(length(tumor) > 0){
  writeFile(tumor, paste0(cancerName, ".rnaseq2.count.tumor.tsv"))
}

totalCounts<-apply(exp, 2, sum)

#per kilobases
genelength<-(posframe$s2 - posframe$s1 + 1) / 1000 

#by column
fpkm<-exp / totalCounts

#by row
fpkm <-t(t(fpkm) / genelength)
fpkm <- fpkm * 1000000

writeFile(fpkm, file=paste0(cancerName, ".rnaseq2.fpkm.tsv"))

normal<-fpkm[,grepl("11$", colnames(fpkm))]
if(is.matrix(normal) && ncol(normal) > 0){
  writeFile(normal, paste0(cancerName, ".rnaseq2.fpkm.normal.tsv"))
}

tumor<-fpkm[,grepl("01$", colnames(fpkm))]
if(length(tumor) > 0){
  writeFile(tumor, paste0(cancerName, ".rnaseq2.fpkm.tumor.tsv"))
}
