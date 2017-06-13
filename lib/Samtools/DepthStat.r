options(bitmapType='cairo')
#resultDir="E:/sqh/Dropbox/career/papers/draft/20160411-guoyan-FFPE"
#inputFile="callable_sites.tsv"
#outputFile="callable_sites.png"
#minimumDepth=20

library(ggplot2)
library(cowplot)

setwd(resultDir)
data<-read.table(inputFile, sep="\t", header=T)

png(file=outputFile, width=2000, height=1500, res=300)
g<-ggplot(data, aes(Group, Count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.5)) +
  xlab("") +
  ylab(paste0("Callable Sites (Depth >= ", minimumDepth, ")")) +
  theme_cowplot()
print(g)
dev.off()
