
filelist<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactor=F)
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
colnames(all)<-filelist$V2
write.csv(all, outFile)
