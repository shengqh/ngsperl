library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)
countFile = args[1]
outputFile = args[2]

#countFile = "/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79.tsv"
#outputFile = "/scratch/cqs/shengq2/vickers/20180410_smallRNA_3018-KCV-77_78_79_mouse_v3/nonhost_genome/bowtie1_bacteria_group1_pm_mappedreads_host_mismatch_table/result/KCV_3018_77_78_79.tsv.pdf"

options(bitmapType='cairo')

counts <- read.table(countFile, sep="\t", header=T, row.names=1)
datForFigure<-t(t(counts)/colSums(counts,na.rm=T))

categoryFigure<-melt(datForFigure)
colnames(categoryFigure) <- c("Category", "Sample", "Count")
categoryFigure$Sample<-factor(categoryFigure$Sample,levels=colnames(counts))

textSize<-9
width<-max(7, ceiling(sqrt(ncol(counts))) * 1.5)

pdf(file=outputFile, width=width, height=width)
p = ggplot(categoryFigure, aes_string(x=factor(1), y="Count", fill="Category")) +
  geom_bar(width=1, stat='identity', color='black') +
  guides(fill=guide_legend(keywidth = 1.5, keyheight = 1.5,override.aes=list(colour=NA))) + # removes black borders from legend
  coord_polar(theta='y') +
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank()) +
  theme(panel.background = element_rect(fill = "white"))+
  theme(legend.text=element_text(size=textSize),
        legend.title= element_text(size=textSize),
        strip.text.x = element_text(size = textSize),
        strip.text.y = element_text(size = textSize)) +
        facet_wrap(~Sample)
print(p)
dev.off()
