rm(list=ls()) 
outFile='P11943_cutrun_hg38'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile6='fileList6.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20240802_11943_cutrun_hg38/bowtie2_summary/result')

### Parameter setting end ###

source("AlignmentUtils.r")
options(bitmapType='cairo')

library(reshape2)
library(ggplot2)
library(stringr)

if (parSampleFile2 != ""){
  draw_chromosome_count(parSampleFile2, outFile)
}

if (parSampleFile1 != ""){
  filelist = read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

  final=NULL
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]
    subdata = read.table(filelocation, sep=":", header=F, strip.white=T, fill=T, stringsAsFactors = F, comment='*')
    totalReads = gsub("\\s.+", "", subdata[1,1])

    subdata = read.table(filelocation, sep=":", header=F, strip.white=T, fill=T, stringsAsFactors = F, comment='*')
    isPairedEnd <- any(grepl("concordantly", subdata$V1))
    totalReads = gsub("\\s.+", "", subdata[1,1])

    if (isPairedEnd) {
      unmappedReads <- subdata$V1[grepl("aligned concordantly 0 times$", subdata$V1)]
      unmappedReads<-gsub("\\s.+", "", unmappedReads)

      uniqueReads <- subdata$V1[grepl("aligned concordantly exactly 1 time$", subdata$V1)]
      uniqueReads<-gsub("\\s.+", "", uniqueReads)

      multiMappedReads <- subdata$V1[grepl("aligned concordantly >1 times$", subdata$V1)]
      multiMappedReads<-gsub("\\s.+", "", multiMappedReads)
    }else{
      unmappedReads <- subdata$V1[grepl("aligned 0 times$", subdata$V1)]
      unmappedReads<-gsub("\\s.+", "", unmappedReads)

      uniqueReads <- subdata$V1[grepl("aligned exactly 1 time$", subdata$V1)]
      uniqueReads<-gsub("\\s.+", "", uniqueReads)

      multiMappedReads <- subdata$V1[grepl("aligned >1 times$", subdata$V1)]
      multiMappedReads<-gsub("\\s.+", "", multiMappedReads)
    }

    if(is.null(final)){
      final = data.frame("TotalReads"=totalReads, "UnmappedReads"=unmappedReads, "UniqueReads"=uniqueReads, "MultiMappedReads"=multiMappedReads, stringsAsFactors = F)
    }else{
      final <-rbind(final, c(totalReads, unmappedReads, uniqueReads, multiMappedReads))
    }
  }
  rownames(final)<-filelist$V2
  write.csv(file=paste0(outFile, ".reads.csv"), final)

  reads=final[,c("UniqueReads", "UnmappedReads", "MultiMappedReads")]
  meltreads=reshape2::melt(t(reads))
  colnames(meltreads)<-c("Read", "Sample", "Count")
  meltreads$Count<-as.numeric(as.character(meltreads$Count))
  meltreads$Read<-factor(meltreads$Read, levels=c("MultiMappedReads", "UnmappedReads", "UniqueReads"))

  sample_height=get_longest_text_width(as.character(meltreads$Sample), "", 11, "inches", 300)

  g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
    geom_bar(stat="identity", width=0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold"),
          axis.text.y = element_text(size=11),
          axis.title.y = element_text(size=11, face="bold"),
          axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size=11, face="bold")) +
    ylab("Reads")
  width=max(2500, 50 * nrow(meltreads)) / 300
  ggsave(file=paste0(outFile, ".reads.png"), g, height=sample_height + 3, width=width, dpi=300, units="in", bg="white")
}
