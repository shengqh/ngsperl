library(SummarizedExperiment)

counts_to_fpkm <- function(counts, featureLength) {
  stopifnot(length(featureLength) == nrow(counts))
  
  # Process one column at a time.
  result <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    countI<-counts[,i]
    totalCount<-sum(countI)
    res<-(10^9)*countI/(featureLength*totalCount)
    return(res)
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(result) <- colnames(counts)
  rownames(result) <- rownames(counts)
  return(result)
}

counts_to_tpm <- function(counts, featureLength) {
  stopifnot(length(featureLength) == nrow(counts))
  
  logFeatureLength<-log(featureLength)
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - logFeatureLength[i]
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

DEBUG=0

args = commandArgs(trailingOnly=TRUE)
cancerName = args[1]
outputFolder = args[2]

if(DEBUG){
  cancerName = "COAD"
  outputFolder = "/scratch/cqs/shengq2/guoyan/prepareRnaseq/result/COAD"
}

setwd(outputFolder)

project = paste0("TCGA-", cancerName)
countFile<-paste0(cancerName, ".rnaseq2.count.rdata")

if(!file.exists(countFile)){
  library(TCGAbiolinks)
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
pos<-data.frame(pos)
pos<-pos[,c("ensembl_gene_id", "seqnames", "start", "end", "gene_id")]

exp<-assay(rnaseq)
rownames(exp)<-pos$ensembl_gene_id

isdup<-duplicated(pos$ensembl_gene_id)
pos<-pos[!isdup,]
exp<-exp[!isdup,]

#sort pos by chromosome and start
pos$length<-pos$end-pos$start+1
pos$seqnames<-substring(pos$seqnames, 4)
pos$seqnamesX<-grepl('[XYZ]', pos$seqnames)
pos$seqnamesN<-as.numeric(pos$seqnames)
pos<-pos[with(pos, order(seqnamesX, seqnamesN, seqnames, start, end)), c("ensembl_gene_id", "seqnames",    "start" ,     "end"      ,"length","gene_id")]

#save position file
genepos<-pos[,c("ensembl_gene_id", "seqnames",    "start" ,     "end" )]
colnames(genepos)<-c("geneid", "chr", "s1", "s2")
write.table(genepos, file=paste0(cancerName, ".rnaseq2.gene.pos"), sep="\t", row.names=F, quote=F)

exp<-exp[as.character(pos$ensembl_gene_id),]
colnames(exp)<-substring(colnames(exp), 1, 15)
colnames(pos)<-c("Feature", "Feature_chr", "Feature_start", "Feature_end", "Feature_length", "Feature_gene_name")

fpkm<-counts_to_fpkm(exp, pos$Feature_length)
tpm<-counts_to_tpm(exp, pos$Feature_length)

normalNames<-colnames(exp)[grepl("11$", colnames(exp))]
tumorNames<-colnames(exp)[grepl("01$", colnames(exp))]

names(normalNames)<-gsub("-11", "", normalNames)
names(tumorNames)<-gsub("-01", "", tumorNames)

pairShortNames<-intersect(names(normalNames), names(tumorNames))
normalPairNames<-normalNames[pairShortNames]
tumorPairNames<-tumorNames[pairShortNames]

exportNames<-list()
exportNames[["normal"]]<-normalNames
exportNames[["tumor"]]<-tumorNames
exportNames[["normal_paired"]]<-normalPairNames
exportNames[["tumor_paired"]]<-tumorPairNames

dataMap<-list()
dataMap[["count"]]<-exp
dataMap[["fpkm"]]<-fpkm
dataMap[["tpm"]]<-tpm

for(dataName in names(dataMap)){
  data<-dataMap[[dataName]]  
  write.table(cbind(pos, data), file=paste0(cancerName, ".rnaseq2.", dataName, ".tsv"), sep="\t", row.names=F, quote=F)
  
  for(exportName in names(exportNames)){
    fileNames<-exportNames[[exportName]]
    subdata<-data[,fileNames]
    if(is.matrix(subdata) && ncol(subdata) > 0){
      write.table(cbind(pos, subdata), file=paste0(cancerName, ".rnaseq2.", dataName, "." , exportName, ".tsv"), sep="\t", row.names=F, quote=F)
    }
  }
}
