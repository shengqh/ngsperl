source("AlignmentUtils.r")

options(bitmapType='cairo')

library(reshape2)
library(ggplot2)

if (parSampleFile3 != '') {
  draw_chromosome_count(parSampleFile3, outFile)
}

if (parSampleFile1 != '') {
  outputPrefix = paste0(outFile, ".STARSummary")
  filelist = read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
  
  final=NULL
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]
    subdata = read.table(filelocation, sep="|", header=F, strip.white=T, fill=T, stringsAsFactors = F, row.names=1)
    colnames(subdata)= (filename)
    if(is.null(final)){
      final = subdata
    }else{
      final = cbind(final, subdata, stringsAsFactors=F)
    }
  }
  
  write.csv(file=paste0(outputPrefix, ".details.csv"), final)
  
  reads=final[c("Number of input reads", "Uniquely mapped reads number", "Number of reads mapped to multiple loci", "Number of reads mapped to too many loci"),]
  rownames(reads)=c("Total", "Unique", "Multiple1", "Multiple2")
  reads[] <- lapply(reads, function(x) {
    as.numeric(as.character(x))
  })
  treads=data.frame(t(reads))
  write.csv(file=paste0(outputPrefix, ".csv"), treads)
  
  treads$Multiple=treads$Multiple1+treads$Multiple2
  treads$Unmapped=treads$Total-treads$Unique-treads$Multiple
  treads=treads[,c(2,5,6)]
  treads$Sample=rownames(treads)
  
  meltreads=melt(treads, id="Sample", variable.name="Read", value.name="Count")
  
  width=max(2000, 50 * nrow(treads))
  png(file=paste0(outputPrefix, ".csv.png"), height=1500, width=width, res=300)
  g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) +
    geom_bar(stat="identity", width=0.5) +
    theme_classic() + xlab("") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
          axis.text.y = element_text(size=11, face="bold"))
  print(g)
  dev.off()
}

if (parSampleFile2 != ''){
  outputPrefix = paste0(outFile, ".FeatureCountSummary")
  filelist = read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
  
  final=NULL
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]
    subdata = read.table(filelocation, sep="\t", header=T, strip.white=T, fill=T, stringsAsFactors = F, row.names=1)
    colnames(subdata)= (filename)
    if(is.null(final)){
      final = subdata
    }else{
      final = cbind(final, subdata)
    }
  }
  
  rowSumCount<-apply(final, 1, sum)
  final<-final[rowSumCount != 0,]
  write.csv(file=paste0(outputPrefix, ".details.csv"), final)
  
  final<-final[rownames(final) != 'Unassigned_MultiMapping',]
  
  treads=data.frame(t(data.matrix(final)))
  treads<-treads[,c(which(!(colnames(treads) == 'Assigned')), which(colnames(treads)=='Assigned'))]
  treads$Total<-apply(treads, 1, sum)
  treads$Percent_Assigned<-treads$Assigned / treads$Total
  
  write.csv(file=paste0(outputPrefix, ".csv"), treads)
  
  tfinal<-data.frame(t(final))
  tfinal$Sample<-factor(rownames(tfinal))
  
  mfinal<-melt(tfinal, id="Sample", variable.name="Read", value.name="Count")
  mfinal$Read<-factor(mfinal$Read, levels=sort(as.character(unique(mfinal$Read))))
  
  width=max(2000, 70 * ncol(final))
  png(file=paste0(outputPrefix, ".csv.png"), height=1500, width=width, res=300)
  
  g<-ggplot(mfinal, aes(x=Sample, y=Count, fill=Read)) +
    geom_bar(stat="identity", width=0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
          axis.text.y = element_text(size=11, face="bold"))
  print(g)
  
  dev.off()
}
