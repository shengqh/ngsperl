rm(list=ls()) 
outFile='AHA_Obesity_EV'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile5='fileList5.txt'
parSampleFile6='fileList6.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20240508_AHA_Emeli_obesity_EV/20240508_rnaseq_hg38/star_featurecount_summary/result')

### Parameter setting end ###

source("AlignmentUtils.r")
options(bitmapType='cairo')

library(reshape2)
library(ggplot2)

if (parSampleFile3 != '') {
  draw_chromosome_count(parSampleFile3, outFile)
  if(exists('parSampleFile6')){
    opt = read.table(parSampleFile6, sep="\t", header=F, stringsAsFactors = F)
    remove_chrM_genes = opt$V1[1] == "1"
    if(remove_chrM_genes){
      draw_chromosome_count(parSampleFile3, paste0(outFile, ".no_chrM"), remove_chrM_genes=TRUE)
    }
  }
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

  fontsize_inch = GeomLabel$default_aes$size * 0.0393701
  
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
  
  g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) +
    geom_bar(stat="identity", width=0.5) +
    theme_classic() + xlab("") +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
          axis.text.y = element_text(size=11, face="bold"))

  width=max(7, fontsize_inch * nrow(treads) + 1)
  ggsave(file=paste0(outputPrefix, ".csv.png"), height=5, width=width, dpi=300, bg="white", limitsize = FALSE)
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
  
  g<-ggplot(mfinal, aes(x=Sample, y=Count, fill=Read)) +
    geom_bar(stat="identity", width=0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
          axis.text.y = element_text(size=11, face="bold"))

  width=max(7, fontsize_inch * ncol(final) + 1)
  ggsave(file=paste0(outputPrefix, ".csv.png"), height=5, width=width, dpi=300, bg="white", limitsize = FALSE)
}

if(exists('parSampleFile5')){
  draw_chromosome_count(parSampleFile5, paste0(outFile, ".gene"))
  draw_gene_count(listFile=parSampleFile5, outFilePrefix=paste0(outFile, ".gene.count"))
  if(exists('parSampleFile6')){
    opt = read.table(parSampleFile6, sep="\t", header=F, stringsAsFactors = F)
    remove_chrM_genes = opt$V1[1] == "1"
    if(remove_chrM_genes){
      draw_chromosome_count(parSampleFile5, paste0(outFile, ".gene.no_chrM"), remove_chrM_genes=TRUE)
    }
  }
}
