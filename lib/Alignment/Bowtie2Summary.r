
options(bitmapType='cairo')

library(reshape2)
library(ggplot2)
library(stringr)

if (parSampleFile1 != ""){
  filelist = read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

  final=NULL
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]
    subdata = read.table(filelocation, sep="\t", header=F, stringsAsFactors = F)
    subdata=subdata[,c(1,3)]
    colnames(subdata)<-c("Chrom", "Reads")
    subdata=subdata[str_length(subdata$Chrom) < 6,]
    subdata=subdata[subdata$Chrom != '*',]
    subdata$Sample=filename
    final=rbind(final, subdata )
  }
  write.csv(file=paste0(outFile, ".chromosome.csv"), final, row.names=F)

  chroms=paste0("chr", c(1:22,'X','Y','M'))
  if(!any(final$Chrom %in% chroms)){
    chroms=paste0("", c(1:22,'X','Y','M'))
  }
  chroms=chroms[chroms %in% final$Chrom]
  final$Chrom=factor(final$Chrom, levels=chroms)
  final$NoRead=final$Reads==0

  colors=c("black","red")
  names(colors)=c("FALSE","TRUE")

  height=max(1000, 100 * length(unique(final$Sample)))
  png(file=paste0(outFile, ".chromosome.png"), height=height, width=3000, res=300)
  g<-ggplot(final, aes(x=Chrom, y=Sample)) + 
    geom_point(aes(size=Reads, color=NoRead)) + theme_classic() + 
    scale_color_manual(values=colors) +
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0),
          axis.title = element_blank())
  print(g)
  dev.off()
}

if (parSampleFile2 != ""){
  filelist = read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)

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

    if(is.null(final )){
      final = data.frame("TotalReads"=totalReads, "UnmappedReads"=unmappedReads, "UniqueReads"=uniqueReads, "MultiMappedReads"=multiMappedReads, stringsAsFactors = F)
    }else{
      final <-rbind(final, c(totalReads, unmappedReads, uniqueReads, multiMappedReads))
    }
  }
  rownames(final)<-filelist$V2
  write.csv(file=paste0(outFile, ".reads.csv"), final)

  reads=final[,c("UniqueReads", "UnmappedReads", "MultiMappedReads")]
  meltreads=melt(t(reads))
  colnames(meltreads)<-c("Read", "Sample", "Count")
  meltreads$Count<-as.numeric(as.character(meltreads$Count))

  width=max(2000, 50 * nrow(meltreads))
  png(file=paste0(outFile, ".reads.png"), height=1500, width=width, res=300)
  g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
    geom_bar(stat="identity", width=0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
          axis.text.y = element_text(size=11, face="bold"))
  print(g)
  dev.off()
}

