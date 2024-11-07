rm(list=ls()) 
outFile='cutrun_11425_mm10_tf'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20240429_cutrun_11425_chalab_mm10_tf/fastqc_raw_summary/result')

### Parameter setting end ###

#outputdir<-"/gpfs23/scratch/cqs/shengq2/rolanda_lister/20190403_rnaseq_2499_lister_mouse_placenta_cutadapt/fastqc_raw/result"
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
  option_tb=read.table(parSampleFile2, sep="\t", header=FALSE)
  myoptions=split(option_tb$V1, option_tb$V2)
  prefix = paste0(myoptions$task_name, ".FastQC")
}else{
  prefix = args[1]
  rmdfile = args[2]
}

summaryfile<-paste0(prefix, ".summary.tsv")
readfile<-paste0(prefix, ".reads.tsv")
baseQualityFile<-paste0(prefix, ".baseQuality.tsv")
sequenceGCFile<-paste0(prefix, ".sequenceGC.tsv")
adapterFile<-paste0(prefix, ".adapter.tsv")
#outputdir<-"H:/shengquanhu/projects/rolanda_lister/20190403_rnaseq_2499_lister_mouse_placenta_cutadapt/temp"
#summaryfile<-"rnaseq_Diabetic_E175.summary.txt"
#readfile<-"rnaseq_Diabetic_E175.reads.txt"
#baseQualityFile<-"rnaseq_Diabetic_E175.baseQuality.txt"
#sequenceGCFile<-"rnaseq_Diabetic_E175.sequenceGC.txt"
#adapterFile<-"rnaseq_Diabetic_E175.adapter.txt"

#setwd(outputdir)

library(ggplot2)
library(reshape2)
library(dplyr)

source("countTableVisFunctions.R")

#summary
fp<-read.table(summaryfile, header=T, sep="\t", stringsAsFactor=F)
fp$QCResult<-factor(fp$QCResult, levels=c("PASS","WARN","FAIL"))
fp$Sample<-factor(fp$Sample, levels=sort(unique(fp$Sample)))
fp$File<-factor(fp$File, levels=sort(unique(fp$File)))
fp$Category<-factor(fp$Category, levels=sort(unique(fp$Category), decreasing=T))

file_height=get_longest_text_width(as.character(fp$File), "", 11, "inches", 300) 
sample_height=get_longest_text_width(as.character(fp$Sample), "", 11, "inches", 300)

g<-ggplot(fp, aes(File, Category))+
  geom_tile(data=fp, aes(fill=QCResult), color="white") +
  scale_fill_manual(values=c("light green", "skyblue", "red")) +
  theme(axis.text.x = element_text(angle=90, vjust=1, size=11, hjust=1, face="bold"),
        axis.text.y = element_text(size=11, face="bold")) +
  xlab("") + ylab("") +
  coord_equal()

width=min(max(2500, 60 * length(unique(fp$File))) + 600, 10000) / 300
ggsave(file=paste0(summaryfile, ".png"), g, height=file_height + 5, width=width, dpi=300, units="in", bg="white")

#reads
fp<-read.table(readfile, header=T, sep="\t") |>
  dplyr::select(Sample, Reads) |>
  dplyr::distinct()

g<-ggplot(fp, aes(x=Sample, y=Reads))+ 
  geom_bar(stat="identity", width=.5)+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold"),
        axis.text.y = element_text(size=11),
        axis.title.y = element_text(size=11, face="bold"),
        axis.title.x = element_blank())

width=min(max(2500, 60 * nrow(fp)), 10000) / 300
ggsave(file=paste0(readfile, ".png"), g, height=sample_height + 3, width=width, dpi=300, units="in", bg="white")

#base quality
fp<-read.table(baseQualityFile, header=T, sep="\t")
fp$newBase = unlist(lapply(fp$Base, function(x){
  if (grepl("-", x)){
    xx=unlist(strsplit(as.character(x), '-'))
    medxx=mean(as.numeric(xx))
    return(medxx)
  }else{
    return(as.numeric(as.character(x)))
  }
}))

xmax=max(fp$newBase)
g<-ggplot(fp, aes(x=newBase, y=Mean)) +
  geom_rect(data=NULL,aes(xmin=0,xmax=xmax, ymin=0,ymax=20),
            fill="darksalmon") +
  geom_rect(data=NULL,aes(xmin=0,xmax=xmax,ymin=20,ymax=28),
            fill="cornsilk") +
  geom_rect(data=NULL,aes(xmin=0,xmax=xmax,ymin=28,ymax=40),
            fill="darkseagreen1") +
  geom_line(aes(group=File), color="limegreen") + xlab("Position in read (bp)") + ylab("Phred score") +
  scale_x_continuous(breaks=seq(0, xmax, 10)) +
  scale_y_continuous(breaks=seq(0, 40, 5)) +
  theme_classic() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = NA),
        panel.grid.major.y = element_line( linewidth=.1, color="gray" ),
        panel.ontop = TRUE)

ggsave(file=paste0(baseQualityFile, ".png"), g, height=1000, width=2500, dpi=300, units="px", bg="white")

#sequence GC
fp<-read.table(sequenceGCFile, header=T, sep="\t")
fp<-group_by(fp, File) %>% mutate(Percent = Count * 100 /sum(Count))
g<-ggplot(fp, aes(x=GC.Content, y=Percent, color=File, group=File)) +
  geom_line() +
  xlab("% GC") + ylab("Percentage of Reads") +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line( linewidth=.1, color="gray" ) )

ggsave(file=paste0(sequenceGCFile, ".png"), g, height=1000, width=2500, dpi=300, units="px", bg="white")

#adapter
fp<-read.table(adapterFile, header=T, sep="\t", quote="", check.names = F)
fpMax<-apply(fp[,c(4:ncol(fp))], 2, max)
fpMaxName = names(fpMax[which(fpMax == max(fpMax))])[1]
colnames(fp)[which(colnames(fp) == fpMaxName)]<-"Adapter"

fp$newPosition = unlist(lapply(fp$Position, function(x){
  if (grepl("-", x)){
    xx=unlist(strsplit(as.character(x), '-'))
    medxx=mean(as.numeric(xx))
    return(medxx)
  }else{
    return(as.numeric(as.character(x)))
  }
}))

ylimmax<-max(50, max(fp$Adapter))
g<-ggplot(fp, aes(x=newPosition, y=Adapter, color=File, group=File)) +
  geom_line() +
  xlab("Position") + ylab("Percentage") + ggtitle(fpMaxName) +
  scale_x_continuous(breaks=seq(0, max(fp$newPosition), 10)) + ylim(0, ylimmax) +
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.major.y = element_line( linewidth=.1, color="gray" ) )

ggsave(file=paste0(adapterFile, ".png"), height=1000, width=2500, dpi=300, units="px", bg="white")

if(length(args) > 0){
  library(knitr)
  rmarkdown::render(rmdfile)
}
