#setwd("/scratch/cqs/shengq1/dnaseq/2110/cnmops/result")
#SampleNames <- c(
#"2110_JP_01"
#,"2110_JP_02"
#,"2110_JP_03"
#,"2110_JP_04"
#,"2110_JP_05"
#,"2110_JP_06"
#,"2110_JP_07"
#,"2110_JP_08"
#,"2110_JP_09"
#,"2110_JP_10"
#,"2110_JP_11"
#,"2110_JP_12"
#,"2110_JP_13"
#,"2110_JP_14"
#,"2110_JP_15"
#,"2110_JP_16"
#)
#BAMFiles <- c(
#"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-1_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-2_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-3_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-4_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-5_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-6_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-7_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-8_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-9_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-10_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-11_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-12_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-13_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-14_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-15_realigned_recal_rmdup.sorted.bam"
#,"/scratch/cqs/shengq1/dnaseq/2110/bwa_markdup/2110-JP-16_realigned_recal_rmdup.sorted.bam"
#)
#hasbed<-1
#bedfile<-"/scratch/cqs/lij17/cnv/SureSelect_XT_Human_All_Exon_V4_withoutchr_withoutY_lite.bed"
#prefix<-"2110"
#callfile<-"2110.call"
#pairmode<-"paired"
#parallel<-8

library(cn.mops)
resfile<-paste0(prefix, "_resCNMOPS_exomecn.mops.Rdata")
if(hasbed){
	segfile<-paste0(prefix, "_x_getSegmentReadCountsFromBAM.Rdata")
	segments <- read.table(bedfile, sep="\t", as.is=TRUE, header=T)
	gr <- GRanges(segments[,1], IRanges(segments[,2],segments[,3]), gene=segments[,4])
	x <- getSegmentReadCountsFromBAM(BAMFiles, GR=gr, sampleNames=SampleNames, mode=pairmode, parallel=parallel)
	save(x, file=segfile)
	resCNMOPS <- exomecn.mops(x)
}else{
	countfile<-paste0(prefix, "_x_getReadCountsFromBAM.Rdata")
	x <- getReadCountsFromBAM(BAMFiles, sampleNames=SampleNames, mode=pairmode) 
	save(x, file=countfile) 
	resCNMOPS <- cn.mops(x) 
}
save(resCNMOPS, file=resfile)
cnvs<-resCNMOPS@cnvs
d<-cbind(as.character(cnvs@elementMetadata@listData$sampleName),
         as.character(cnvs@seqnames),
         as.character(cnvs@ranges@start),
         as.character(as.numeric(cnvs@ranges@start) + as.numeric(cnvs@ranges@width) - 1),
         as.character(cnvs@ranges@width),
         as.character(cnvs@elementMetadata@listData$CN),
         as.character(cnvs@elementMetadata@listData$median),
         as.character(cnvs@elementMetadata@listData$mean))
colnames(d)<-c("sample","chr","start","end", "length","type","median","mean")
d<-d[order(d[,"sample"], as.numeric(d[,"chr"]), as.numeric(d[,"start"])),]
d[,"chr"]<-paste0("chr",d[,"chr"])
d[,"type"]<-apply(d,1,function(x){
  if(as.numeric(x["median"]) < 0){
    return ("DELETION")
  }else{
    return ("DUPLICATION")
  }
})
write.table(d, file=callfile,sep="\t",col.names=T,row.names=F,quote=F)
