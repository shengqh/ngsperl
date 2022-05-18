
options(bitmapType='cairo')

library(testit)
library(ChIPpeakAnno)
genes <- read.table(parFile1, sep="\t", header=F)
dup<-genes[duplicated(genes$V4),]
dup$V4<-paste0(dup$V4,":",dup$V2)
genes[rownames(dup), "V4"] = dup$V4

assert(length(duplicated(genes$V4)) == 0)
 
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
  
  write.table(df, file=paste0(file, ".nearest_gene.txt"), quote=F, row.names=F, sep="\t")

  res=rbind(res, data.frame("File"=paste0(name, ".nearest_gene.txt"), "Sample"=name))
}

write.csv(res, paste0(outFile, ".csv"), row.names=F)
