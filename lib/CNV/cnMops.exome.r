rm(list=ls()) 
outFile='P11796'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/cqs/references/exomeseq/Twist/Twist_Mouse_Exome_Target_Rev1_7APR20.main_chr.bed'
parFile2=''
parFile3=''

setwd('/nobackup/h_vangard_1/breast_cancer_spore/20240724_balko_11796_exome_mm10/bwa_g4_refine_cnMOPS/result')

### Parameter setting end ###

options(bitmapType='cairo')

library(data.table)
library(GenomicRanges)
library(cn.mops)
library(DNAcopy)
library(dplyr)

bam_map=fread(parSampleFile1, header=FALSE, data.table=FALSE)
sample_names=bam_map$V2
bam_files=bam_map$V1

bedfile<-parFile1

myoptions_tbl=fread(parSampleFile2, header=FALSE)
myoptions=split(myoptions_tbl$V1, myoptions_tbl$V2)
parallel<-as.numeric(myoptions$thread)

refnames=fread(parSampleFile3, header=FALSE)$V1
if(all(is.na(refnames))){
  refnames=c()
}

segments <- read.table(bedfile, sep="\t", as.is=TRUE, header=T)
gr <- GRanges(segments[,1], IRanges(segments[,2]-30, segments[,3]+30))
gr <- reduce(gr)
sort(gr, ignore.strand=TRUE)

if(length(refnames) > 0){
  insample<-sample_names %in% refnames
  REFNames<-sample_names[insample]
  REFFiles<-bam_files[insample]
  SAMNames<-sample_names[!insample]
  SAMFiles<-bam_files[!insample]

  segfile<-paste0(outFile, ".getSegmentReadCountsFromBAM_ref.Rdata")
  if(file.exists(segfile)){
    load(segfile)
  }else{
    refdata <- getSegmentReadCountsFromBAM(REFFiles, GR=gr, sampleNames=REFNames, parallel=parallel)
    samdata <- getSegmentReadCountsFromBAM(SAMFiles, GR=gr, sampleNames=SAMNames, parallel=parallel)
    save(refdata, samdata, file=segfile)
  }

  resCNMOPS<-referencecn.mops(cases=samdata, 
                              controls=refdata)
}else{
  segfile<-paste0(outFile, ".getSegmentReadCountsFromBAM.Rdata")
  if(file.exists(segfile)){
    load(segfile)
  }else{
    x <- getSegmentReadCountsFromBAM(bam_files, GR=gr, sampleNames=sample_names, parallel=parallel)
    save(x, file=segfile)
  }
  resCNMOPS<-exomecn.mops(x)
}

resCNMOPS<-calcIntegerCopyNumbers(resCNMOPS)
resfile<-paste0(outFile, ".resCNMOPS.rds")
saveRDS(resCNMOPS, file=resfile)

#individual sample level
d<-as.data.frame(cnvs(resCNMOPS)) |>
  dplyr::mutate(type=ifelse(median < 0, "DELETION", "DUPLICATION"),
                chrNumStr=gsub("Y", "24", gsub("X", "23", gsub("chr","",seqnames))),
                chrNum=as.numeric(chrNumStr)) |>
  dplyr::arrange(sampleName, chrNum, start) |>
  dplyr::select(-chrNumStr)
write.table(d |> dplyr::select(-chrNum), file=paste0(outFile, ".cnMops.tsv"),sep="\t",col.names=T,row.names=F,quote=F)

locus<-d |> 
  dplyr::mutate(name=paste0(seqnames, ".", start, ".", end, ".", CN, ".", sampleName)) |>
  dplyr::arrange(seqnames, start, sampleName) |>
  dplyr::select(seqnames, start, end, name)
write.table(locus, file=paste0(outFile, ".cnMops.bed"), sep="\t", col.names=F, row.names=F,quote=F)

