#https://bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#transformations-from-the-raw-scale

library(tools)

args = commandArgs(trailingOnly=TRUE)

print(paste0("length of args = ", length(args)))

if(length(args) == 0){
  inputFile = "/data/cqs/justin_turner_data/20210520_CRS/CRS_Swabs_Hisat2Counts.csv"
  outputFile = "/data/cqs/justin_turner_data/20210520_CRS/CRS_Swabs_Hisat2Counts.gene.csv"
  gtfMapFile = "/data/cqs/references/gencode/GRCh38.p13/gencode.v38.annotation.gtf.map"
}else{
  inputFile = args[1]
  outputFile = args[2]
  gtfMapFile = args[3]
}
print(paste0("inputFile=", inputFile))
print(paste0("outputFile=", outputFile))
print(paste0("gtfMapFile=", gtfMapFile))

if(file_ext(inputFile) == "csv"){
  dat<-read.csv(inputFile, row.names=1)
}else{
  dat<-read.table(inputFile, row.names=1, sep="\t", header=T)
}

ignoreColumns = c("geneName", "geneDescription")
dat<-dat[,!(colnames(dat) %in% ignoreColumns)]
dat<-dat[rowSums(dat) >0, ]
nsample=ncol(dat)

gtfMap<-read.table(gtfMapFile, sep="\t", header=T)
gtfMap<-gtfMap[gtfMap$gene_biotype == "protein_coding",]
gtfMap$id<-gsub('\\..+',"",gtfMap$gene_id)
gtfMap<-gtfMap[gtfMap$id %in% rownames(dat), c("id", "gene_name")]
gtfMap<-gtfMap[!duplicated(gtfMap),]
gtfNameMap<-split(gtfMap$gene_name, gtfMap$id)

dat<-dat[rownames(dat) %in% gtfMap$id,]
dat$gene_name<-unlist(gtfNameMap[rownames(dat)])

library(plyr)

dup_dat<-dat[duplicated(dat$gene_name),]
all_dup_dat<-dat[dat$gene_name %in% dup_dat$gene_name,]

unique_dat<-dat[!(dat$gene_name %in% dup_dat$gene_name),]

rm_dup_dat=data.frame()
gene=dup_dat$gene_name[1]
for(gene in unique(dup_dat$gene_name)){
  sub_dat=all_dup_dat[all_dup_dat$gene_name == gene,]
  merged_dat=colSums(sub_dat[,c(1:nsample)])
  rm_dup_dat<-rbind(rm_dup_dat, merged_dat)
}
rm_dup_dat$gene_name=unique(dup_dat$gene_name)
colnames(rm_dup_dat)<-colnames(dat)

all_dat<-rbind(unique_dat, rm_dup_dat)
rownames(all_dat)<-all_dat$gene_name
all_dat<-all_dat[,c(1:nsample)]
write.csv(all_dat, file=outputFile, quote=F)
