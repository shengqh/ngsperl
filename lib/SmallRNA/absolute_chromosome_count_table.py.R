library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
countFile = args[1]
outputFile = args[2]

#setwd("c:/temp")
#countFile = "6369_ES.count.txt"
#outputFile = "6369_ES.count.txt.png"

options(bitmapType='cairo')

counts <- read.table(countFile, sep="\t", header=T, row.names=1)
counts<-counts[grepl("_read", rownames(counts)),]
rownames(counts)<-gsub("_read", "", rownames(counts))

nsample=ncol(counts)

totalreads<-unlist(counts['total',])

counts$Genome<-rownames(counts)

mcounts<-melt(counts, id="Genome")
mcounts<-mcounts[mcounts$Genome != "total",]
mcounts$Perc<-mcounts$value / totalreads[mcounts$variable]

colnames(mcounts)<-c("Genome", "Sample", "Count", "Perc")

nwidth=ceiling(sqrt(nsample))
nheight=ceiling(nsample/nwidth)

png(outputFile, width=max(500 * nwidth, 1000) + 200, height=max(500 * nheight), res=300)
g<-ggplot(mcounts, aes(x=Genome, y=Perc, fill=Genome)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~Sample) + 
  ylab("Percentage") + 
  theme_bw() + 
  theme(strip.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(g)
dev.off()

