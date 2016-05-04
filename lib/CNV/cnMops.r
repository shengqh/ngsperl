#setwd("/scratch/cqs/shengq1/dnaseq/2110/cnmops/result")
#sample_names <- c(
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
#bam_files <- c(
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
#refnames<-c()

if(hasbed){
  segments <- read.table(bedfile, sep="\t", as.is=TRUE, header=T)
  gr <- GRanges(segments[,1], IRanges(segments[,2], segments[,3]))
  gr <- reduce(gr)
  sort(gr, ignore.strand=TRUE)
}    

library(cn.mops)
library(DNAcopy)
resfile<-paste0(prefix, "_resCNMOPS.cnmops.Rdata")

if(length(refnames) > 0){
  insample<-sample_names %in% refnames
  REFNames<-sample_names[insample]
  REFFiles<-bam_files[insample]
  SAMNames<-sample_names[!insample]
  SAMFiles<-bam_files[!insample]
  if(hasbed){
    segfile<-paste0(prefix, "_getSegmentReadCountsFromBAM_ref.Rdata")
    if(file.exists(segfile)){
      load(segfile)
    }else{
      refdata <- getSegmentReadCountsFromBAM(REFFiles, GR=gr, sampleNames=REFNames, mode=pairmode, parallel=parallel)
      samdata <- getSegmentReadCountsFromBAM(SAMFiles, GR=gr, sampleNames=SAMNames, mode=pairmode, parallel=parallel)
      save(refdata, samdata, file=segfile)
    }
  }else{
    countfile<-paste0(prefix, "_getReadCountsFromBAM_ref.Rdata")
    if(file.exists(countfile)){
      load(countfile)
    }else{
      refdata <- getReadCountsFromBAM(REFFiles, sampleNames=REFNames, mode=pairmode, parallel=parallel)
      samdata <- getReadCountsFromBAM(SAMFiles, sampleNames=SAMNames, mode=pairmode, parallel=parallel)
      save(refdata, samdata, file=countfile)
    }
  }
  resCNMOPS<-referencecn.mops(cases=samdata, 
                              controls=refdata, 
                              upperThreshold=0.5, 
                              lowerThreshold=-0.5,
                              segAlgorithm="fast")
  resCNMOPS<-calcIntegerCopyNumbers(resCNMOPS)
}else{
  if(hasbed){
    segfile<-paste0(prefix, "_getSegmentReadCountsFromBAM.Rdata")
    if(file.exists(segfile)){
      load(segfile)
    }else{
      x <- getSegmentReadCountsFromBAM(bam_files, GR=gr, sampleNames=sample_names, mode=pairmode, parallel=parallel)
      save(refx, samx, file=segfile)
    }
    resCNMOPS<-exomecn.mops(x, upperThreshold=0.5, lowerThreshold=-0.5)
    resCNMOPS<-calcIntegerCopyNumbers(resCNMOPS)
  }else{
    countfile<-paste0(prefix, "_getReadCountsFromBAM.Rdata")
    if(file.exists(countfile)){
      load(countfile)
    }else{
      x <- getReadCountsFromBAM(bam_files, sampleNames=sample_names, mode=pairmode)
      save(x, file=countfile)
    }
    resCNMOPS <- cn.mops(x) 
    resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
  }
}

#load(resfile)
save(resCNMOPS, file=resfile)

d<-as.data.frame(cnvs(resCNMOPS))
d[,"type"]<-apply(d,1,function(x){
  if(as.numeric(x["median"]) < 0){
    return ("DELETION")
  }else{
    return ("DUPLICATION")
  }
})

# locus<-data.frame(Title=paste0(d$sampleName, " ~ ", d$seqnames, ":", d$start, "-", d$end, " ~ ", d$type),
#                   Filename=paste0(d$sampleName, "_", d$seqnames, "_", d$start, "_", d$end, ".png"))

d<-d[order(d[,"sampleName"], as.numeric(d[,"seqnames"]), as.numeric(d[,"start"])),]
write.table(d, file=callfile,sep="\t",col.names=T,row.names=F,quote=F)

locus<-d[,c("seqnames", "start", "end")]
locus$name<-paste0(d$seqnames, "_", d$start, "_", d$end, "_", d$CN, "_", d$sampleName)
locus<-locus[order(d$seqnames, d$start),]
write.table(locus, file=paste0(prefix, ".call.bed"), sep="\t", col.names=F, row.names=F,quote=F)

cnvr<- data.frame(seqnames=seqnames(resCNMOPS@cnvr),
                  starts=start(resCNMOPS@cnvr)-1,
                  ends=end(resCNMOPS@cnvr),
                  file=paste(seqnames(resCNMOPS@cnvr), start(resCNMOPS@cnvr)-1, end(resCNMOPS@cnvr), sep="_") )
cnvr<-data.frame(cbind(cnvr, elementMetadata(resCNMOPS@cnvr)))
write.table(file=paste0(prefix, ".cnvr.tsv"), cnvr, sep="\t" ,row.names=F, quote=F)

# dir.create("images", showWarnings = FALSE)
# 
# index<-9
# for(index in c(1:nrow(locus))){
#   png(filename=paste0("images/", locus$Filename[index]), width=2000, height=2000, res=300)
#   plot(resCNMOPS, which=index, toFile=TRUE)
#   grid.text(locus$Title[index], y=0.98)
#   dev.off()
# }
