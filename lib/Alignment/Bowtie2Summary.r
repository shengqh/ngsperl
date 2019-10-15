
options(bitmapType='cairo')

inputFile<-parSampleFile1
outputFile<-outFile

#setwd('/scratch/cqs/shengq1/reiter/20170404_parclip_3495_human/bowtie1/result/')
#inputFile = "/scratch/cqs/shengq1/reiter/20170404_parclip_3495_human/bowtie1/result/fileList1.txt"
#outputFile = "parclip_human"

library(reshape2)
library(ggplot2)

cat("inputFile=", inputFile, "\n")
cat("outputFile=", outputFile, "\n")

filelist = read.table(inputFile, sep="\t", header=F, stringsAsFactors = F)

final=NULL
i=1
for(i in c(1:nrow(filelist))){
  filename = filelist$V2[i]
  filelocation =filelist$V1[i]
  subdata = read.table(filelocation, sep=":", header=F, strip.white=T, fill=T, stringsAsFactors = F, comment='*')
  totalReads = gsub("\\s.+", "", subdata[1,1])
  
  unmappedReads <- subdata$V1[grepl("aligned concordantly 0 times$", subdata$V1)]
  unmappedReads<-gsub("\\s.+", "", unmappedReads)

  uniqueReads <- subdata$V1[grepl("aligned concordantly exactly 1 time$", subdata$V1)]
  uniqueReads<-gsub("\\s.+", "", uniqueReads)
  
  multiMappedReads <- subdata$V1[grepl("aligned concordantly >1 times$", subdata$V1)]
  multiMappedReads<-gsub("\\s.+", "", multiMappedReads)

  if(is.null(final )){
    final = data.frame("TotalReads"=totalReads, "UnmappedReads"=unmappedReads, "UniqueReads"=uniqueReads, "MultiMappedReads"=multiMappedReads, stringsAsFactors = F)
  }else{
    final <-rbind(final, c(totalReads, unmappedReads, uniqueReads, multiMappedReads))
  }
}
rownames(final)<-filelist$V2
write.csv(file=paste0(outputFile, ".csv"), final)

reads=final[,c("UniqueReads", "UnmappedReads", "MultiMappedReads")]
meltreads=melt(t(reads))
colnames(meltreads)<-c("Read", "Sample", "Count")
meltreads$Count<-as.numeric(as.character(meltreads$Count))

width=max(2000, 50 * nrow(meltreads))
png(file=paste0(outputFile, ".png"), height=1500, width=width, res=300)
g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
  geom_bar(stat="identity", width=0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
        axis.text.y = element_text(size=11, face="bold"))
print(g)
dev.off()
  
