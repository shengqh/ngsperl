library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

inputFile = args[1]
outputFile = args[2]

dat<-read.table(inputFile, sep="\t", header=T, stringsAsFactors = F)
dat2<-dat[!(dat$Chromosome %in% c("X", "Y")),]

dmedian<-dat2%>%
  group_by(Chromosome)%>% 
  summarise(Median=median(Count))

medCount<-dmedian$Median
names(medCount)<-dmedian$Chromosome

dat2$SizeFactor<-medCount[dat2$Chromosome] / dat2$Count

smedian<-dat2%>%
  group_by(Sample)%>% 
  summarise(Median=median(SizeFactor))
medSizeFactor<-smedian$Median
names(medSizeFactor)<-smedian$Sample

datx<-dat[dat$Chromosome == "X",]
datx$SizeFactor<-medSizeFactor[datx$Sample]

daty<-dat[dat$Chromosome == "Y",]
daty$SizeFactor<-medSizeFactor[daty$Sample]

finalDat<-rbind(dat2, datx, daty)

finalDat$LogSizeFactor<-log2(finalDat$SizeFactor)

write.table(finalDat, file=outputFile, sep="\t", row.names=F, quote=F)

nsample=ceiling(sqrt(length(unique(finalDat$Sample))))
swidth=max(1000, nsample * 500) 
png(paste0(outputFile, ".png"), width=swidth, height=800, res=300)
g<-ggplot(finalDat, aes(x=Sample,y=SizeFactor)) + geom_violin() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(g)
dev.off()