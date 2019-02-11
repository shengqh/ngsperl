#setwd("")
#covfile<-""
#sample_names<-c("","")

library(ggplot2)
library(reshape2)
library(RColorBrewer)

dat<-read.table(covfile, sep="\t", header=F)
chrom<-dat$V1[1]
dat<-dat[,c(2,4:ncol(dat))]
colnames(dat)<-c("start", sample_names)

csum<-apply(dat[,-1], 1, sum)
csumz<-which(csum>0)
mincsumz<-min(csumz)
maxcsumz<-max(csumz)

cdat<-dat[c(mincsumz:maxcsumz),]
mdat<-melt(cdat, id="start")
colnames(mdat)<-c("start", "Sample", "Reads")
mdat$logCount<-log2(mdat$Reads+1)

makeColors<-function(n,colorNames="Set1") {
  maxN<-brewer.pal.info[colorNames,"maxcolors"]
  if (n<=maxN) {
    colors<-brewer.pal(n, colorNames)
  } else {
    colors<-colorRampPalette(brewer.pal(maxN, colorNames))(n)
  }
  return(colors)
}

sampleColors<-makeColors(ncol(dat)-1)

pdf(file=paste0(covfile, ".pdf"), width=10, height=5)

g<-ggplot(mdat, aes(x=start, y=Reads)) + 
  geom_point(aes(col=Sample), size=0.5) + 
  scale_color_manual(values=sampleColors) +
  xlab(paste("Chromosome", chrom)) +
  ylab("Reads") +
  theme_classic()
print(g)
dev.off()


pdf(file=paste0(covfile, ".log2.pdf"), width=10, height=5)

g<-ggplot(mdat, aes(x=start, y=logCount)) + 
  geom_point(aes(col=Sample), size=0.5) + 
  scale_color_manual(values=sampleColors) +
  xlab(paste("Chromosome", chrom)) +
  ylab("Reads") +
  theme_classic()
print(g)
dev.off()
