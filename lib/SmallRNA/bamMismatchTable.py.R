library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
countFile = args[1]
outputFile = args[2]

#countFile = "/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/bowtie1_nonhost_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79.tsv"
#outputFile = "/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/data_visualization/bowtie1_nonhost_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79.tsv.pdf"

options(bitmapType='cairo')

counts <- read.table(countFile, sep="\t", header=T, row.names=1)

countFigure<-melt(data.matrix(counts))
colnames(countFigure) <- c("Category", "Sample", "Count")
countFigure$Sample<-factor(countFigure$Sample,levels=colnames(counts))

datForFigure<-t(t(counts)/colSums(counts,na.rm=T))
categoryFigure<-melt(datForFigure)
colnames(categoryFigure) <- c("Category", "Sample", "Count")
categoryFigure$Sample<-factor(categoryFigure$Sample,levels=colnames(counts))

countFigure$Type = "Count"
categoryFigure$Type = "Percentage"

allFigure<-rbind(countFigure, categoryFigure)

textSize<-20
xSize<-12
width<-max(8, ceiling(sqrt(ncol(counts))) * 1.5)

pdf(file=outputFile, width=width, height=8, onefile=T)
p = ggplot(allFigure, aes_string(x="Sample", y="Count", fill="Category")) +
  geom_bar(width=1, stat='identity', color='black') +
  facet_wrap(~Type, scales = "free_y", strip.position="left", ncol=1) +
  guides(fill=guide_legend(keywidth = 1.5, keyheight = 1.5,override.aes=list(colour=NA))) + # removes black borders from legend
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size = xSize),
        legend.text=element_text(size=textSize),
        legend.title= element_text(size=textSize),
        legend.position = "top",
        strip.text.x = element_text(size = textSize),
        strip.text.y = element_text(size = textSize),
        strip.placement = "outside",
        strip.background=element_blank()) + 
  ylab("") +
  xlab("")
print(p)

dev.off()
