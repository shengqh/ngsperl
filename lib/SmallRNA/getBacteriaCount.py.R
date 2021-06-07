
library(ggplot2)
library(reshape2)

countFile = parFile1
taskReadFile = parFile2
outputFilePrefix = outFile

#countFile = "/scratch/cqs/kasey_vickers_projects/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76.tsv.summary"
#taskReadFile = "/scratch/cqs/kasey_vickers_projects/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/reads_in_tasks/result/HDL_76.NonParallel.TaskReads.csv"
#outputFilePrefix = "/scratch/cqs/kasey_vickers_projects/20191112_smallRNA_3018-KCV_76_mouse_v4_tRNA_byTiger/data_visualization/bacteria_count/result/HDL_76.tsv.summary"

options(bitmapType='cairo')

counts <- read.table(countFile, sep="\t", header=T, row.names=1, check.names=F)
taskCounts <- data.frame(t(read.csv(taskReadFile, header=T, row.names=1, stringsAsFactor=F, check.names=F)))

totalCounts<-rowSums(taskCounts)
counts$HostSmallRNA<-taskCounts[rownames(counts), "Host.Small.RNA"]

counts$Other<-totalCounts[rownames(counts)] - counts$Count - counts$HostSmallRNA

colnames(counts)<-c("Bacteria", "Host Small RNA", "Other")

countFigure<-melt(data.matrix(counts))
colnames(countFigure) <- c("Sample", "Category",  "Count")
countFigure$Sample<-factor(countFigure$Sample,levels=rownames(counts))
countFigure$Type = "Count"

perc<-data.frame(counts/rowSums(counts,na.rm=T), check.names=F)
percFigure<-melt(data.matrix(perc))
colnames(percFigure) <- c("Sample", "Category",  "Count")
percFigure$Sample<-factor(percFigure$Sample,levels=rownames(counts))
percFigure$Type = "Percentage"

rpm<-data.frame(counts/rowSums(counts,na.rm=T) * 1000000, check.names=F)
write.csv(rpm, file=paste0(outputFilePrefix, ".rpm.csv"))

rpmFigure<-melt(data.matrix(rpm))
colnames(rpmFigure) <- c("Sample", "Category",  "Count")
rpmFigure$Sample<-factor(rpmFigure$Sample,levels=rownames(counts))
rpmFigure$Type = "Reads per million"

allFigure<-rbind(countFigure, rpmFigure)

save(counts, perc, rpm, file=paste0(outputFilePrefix, ".rdata") )

textSize<-20
xSize<-12
width<-max(10, nrow(counts) * 0.3)

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

pdf(file=paste0(outputFilePrefix, ".pdf"), width=width, height=8, onefile=T)
print(p)
dev.off()

png(file=paste0(outputFilePrefix, ".png"), width=max(4000, nrow(counts) * 100), height=3200, res=300)
print(p)
dev.off()
