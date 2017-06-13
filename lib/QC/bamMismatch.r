options(bitmapType='cairo')

library(ggplot2)
library(cowplot)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
inputfile = args[1]
outputfile = args[2]
xlim=as.numeric(args[3])

count = read.table(inputfile, sep="\t", header=T)
percCount = data.frame(count %>% group_by(File) %>%  mutate(Percentage=round(100 * ReadCount / sum(ReadCount), 2)) %>% ungroup)

numberOfSamples=length(unique(count$File))
sqrtOfSample = ceiling(sqrt(numberOfSamples))

if(length(args)>=5){
  height=as.numeric(args[4])
  width=as.numeric(args[5])
}else{
  height=600 * sqrtOfSample
  width=600 * sqrtOfSample
}

png(file=outputfile, height=height, width=width, res=300)
g=ggplot(data=percCount, aes(x=NumberOfMismatch, y=Percentage)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~File) + 
  xlab("Number of mismatch")
		
if(xlim > 0){
  g = g + scale_x_continuous(limits=c(0,xlim), breaks = c(0:xlim))
}

print(g)
dev.off()
