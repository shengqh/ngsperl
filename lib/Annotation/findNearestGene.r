rm(list=ls()) 
outFile='P11057_cutrun_hg38'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/cqs/references/gencode/GRCh38.p13/gencode.v38.annotation.gtf.map.bed'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/shengq2/20240202_11057_cutruntools2_hg38_test/seacr_gene/result')

### Parameter setting end ###

options(bitmapType='cairo')

library(testit)
library(ChIPpeakAnno)

genes <- read.table(parFile1, sep="\t", header=F)
dup<-genes[duplicated(genes$V4),]
dup$V4<-paste0(dup$V4,":",dup$V2)
genes[rownames(dup), "V4"] = dup$V4

assert(sum(duplicated(genes$V4)) == 0)
 
colnames(genes)<-c("seqnames", "start", "end", "name")
gene2<-toGRanges(genes)

files<-read.table(parSampleFile1, sep="\t", header=F)

res=NULL
i=1
for (i in c(1:nrow(files))){
  file = files[i,1]
  name = files[i,2]

  macsOutput <- toGRanges(file, format="BED", skip=1)
  annotated <- annotatePeakInBatch(macsOutput, AnnotationData=gene2)
  
  df<-mcols(annotated)
  df$seqnames<-seqnames(annotated)
  df$peak_start<-start(annotated)
  df$peak_end<-end(annotated)
  
  df<-df[,c("seqnames", "peak_start", "peak_end", "peak", "score", "feature", "start_position", "end_position", "insideFeature", "distancetoFeature", "shortestDistance")]
  colnames(df)<-c("chr", "start", "end", "feature", "score", "gene", "gene_start", "gene_end", "insideGene", "distanceToGene", "shortestDistanceToGene")
  
  gene_file=paste0(file, ".nearest_gene.txt")
  write.table(df, file=gene_file, quote=F, row.names=F, sep="\t")

  res=rbind(res, data.frame("File"=gene_file, "Sample"=name))
}

write.csv(res, paste0(outFile, ".files.csv"), row.names=F)

