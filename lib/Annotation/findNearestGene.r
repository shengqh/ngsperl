setwd("Z:/Shared/Labs/Brown,J/tiger/20170119_atacseq_comparison/bwa_macs2diff/result/HAEC_10perc_vs_5perc")
library(ChIPpeakAnno)

macsOutput <- toGRanges("HAEC_10perc_vs_5perc_c3.0_cond1.bed", format="BED", skip=1)
genes <- toGRanges("H:/shengquanhu/projects/database/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.map.bed", format="BED", skip=1)
annotated <- annotatePeakInBatch(macsOutput, AnnotationData=genes)

df<-mcols(annotated)
df$seqnames<-seqnames(annotated)
df$peak_start<-start(annotated)
df$peak_end<-end(annotated)

df<-df[,c("seqnames", "peak_start", "peak_end", "peak", "score", "feature", "start_position", "end_position", "insideFeature", "distancetoFeature", "shortestDistance")]
colnames(df)<-c("chr", "start", "end", "feature", "score", "gene", "gene_start", "gene_end", "insideGene", "distanceToGene", "shortestDistanceToGene")
write.csv(df, file="test.csv", quote=F, row.names=F)
