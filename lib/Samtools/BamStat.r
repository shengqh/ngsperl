library(dplyr)
library(reshape2)
library(ggplot2)

options(bitmapType='cairo')

if (parSampleFile2 != ""){
  draw_chromosome_count(parSampleFile2, outFile)
}

if (parSampleFile1 != ""){
  filelist<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactor=F)
  if(nrow(filelist) != length(unique(filelist$V2))){
    stop(paste0("Cannot have replicated sample names in ", parSampleFile1, " ."))
  }

  filecounts<-apply(filelist, 1, function(x){
    dat<-read.table(x[1], sep="\t")
    #dat<-read.table(filelist[1,1], sep="\t")
    counts<-sapply(strsplit(as.vector(dat$V1), ' \\+ '), "[", 1)
    names<-sapply(strsplit(as.vector(dat$V1), ' \\+ '), "[", 2)
    names<-sapply(strsplit(names, ' \\('), "[", 1)
    names<-sub("\\S+\\s+", "", names)
    counts<-counts[1:11]
    names(counts)<-names[1:11]
    return(counts)
  })
  all<-rbind(filecounts)
  all<-data.frame(all)
  colnames(all)<-filelist$V2
  write.csv(all, paste0(outFile,".reads.csv"))

  if("properly paired" %in% rownames(all)){
    reads=data.frame(t(all[c("in total", "mapped", "properly paired"),]))
    idat=mutate_all(reads, function(x) as.numeric(as.character(x)))
    colnames(idat)=c("total","mapped","paired")
    idat$unmapped=idat$total-idat$mapped
    idat$unpaired=idat$mapped-idat$paired
    idat=idat[,c("paired","unpaired","unmapped")]
    colnames(idat)=c("Mapped properly paired", "Mapped not properly paired", "Unmapped")
    idat$sample=rownames(idat)
    mdat=melt(idat,id="sample")
    colnames(mdat)=c("Sample", "Reads", "Count")
    png(paste0(outFile, ".reads.png"), width=3000, height=2000, res=300)
    g<-ggplot(mdat, aes(fill=Reads, y=Count, x=Sample)) + 
      geom_bar(position="stack", stat="identity") + theme_bw() + theme(legend.position="top")
    print(g)
    dev.off()
  }
}