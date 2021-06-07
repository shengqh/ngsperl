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
  subdata<-subset(subdata, V2 != "")
  if(nrow(subdata) != 4){
    stop(paste0("The valid data lines are not equals to 4 in file ", filelocation))
  }
  subdata$V2<-sub('\\s+.*', '', subdata$V2)
  colnames(subdata)[2]<-filename
  if(is.null(final)){
    final = subdata[,2,drop=F]
  }else{
    final = cbind(final, subdata[,2,drop=F])
  }
}
row.names(final)<-c("Total", "Mapped", "Failed", "Suppressed")
write.csv(file=paste0(outputFile, ".csv"), final)

reads=final[c("Mapped", "Failed", "Suppressed"),]
meltreads=melt(t(reads))
colnames(meltreads)<-c("Sample", "Read", "Count")
meltreads$Count<-as.numeric(as.character(meltreads$Count))
meltreads$Read<-factor(meltreads$Read, levels=c("Suppressed", "Failed", "Mapped"))

width=max(2000, 50 * nrow(meltreads))
png(file=paste0(outputFile, ".png"), height=1500, width=width, res=300)
g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
  geom_bar(stat="identity", width=0.5) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
        axis.text.y = element_text(size=11, face="bold"))
print(g)
dev.off()
  
