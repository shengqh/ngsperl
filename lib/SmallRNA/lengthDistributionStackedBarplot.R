
source("countTableVisFunctions.R")

library(ggplot2)
library(data.table)
library(reshape2)

fileList1 <- read.table("fileList1.txt", header=FALSE, sep="\t", stringsAsFactors = F, row.names=2)
RNA_class <- rownames(fileList1)[!(rownames(fileList1) %in% c("category"))]

category <- read.csv(fileList1["category",1], header=T, row.names=1)
totalReads<-unlist(category["TotalReads",])

final<-NULL
for(rna in RNA_class){
  print(rna)
  read.count<-fread(fileList1[rna,1],data.table=FALSE,check.names=TRUE)

  if(is.numeric(read.count[,1])){
    read.count.length=read.count[,1]
  }else{
    read.count.length=nchar(read.count[,1])
  }

  if("TopFeature" %in% colnames(read.count)){
    read.count<-read.count[,c(-1,-2,-3)]
  }else{
    read.count<-read.count[,-1]
  }

  read.count.agg<-aggregate(read.count, by=list(read.count.length),sum)
  row.names(read.count.agg)<-read.count.agg[,1]
  read.count.agg<-read.count.agg[,colnames(category)]
  read.count.agg.rpm<-t(t(read.count.agg) / totalReads * 1e+6)
  
  mdata<-reshape2::melt(read.count.agg.rpm)
  colnames(mdata)<-c("Length", "Sample", "RPM")
  mdata$Category<-rna
  if(is.null(final)){
    final<-mdata
  }else{
    final<-rbind(final, mdata)
  }
}

fastq_length<-final[final$Category=="fastq_len",]
final<-final[! (final$Category %in% c("fastq_len", "Genome")),]

RNA_class<-RNA_class[RNA_class != "fastq_len"]
final$Category<-factor(final$Category, levels=(RNA_class))

allcolors<-c("blue","green","red", "brown","purple", "yellow", "black")[1:length(unique(final$Category))]

textTitle<-element_text(face= "bold", color = "black", size=22, hjust=0.5)
text20Bold<-element_text(face= "bold", color = "black", size=20)
text20<-element_text(color = "black", size=20)

fastq_length<-fastq_length[fastq_length$Length >= 10,]

pdf(file=paste0(outFile,".pdf"), onefile=T, width=8, height=7)

for(sample in unique(fastq_length$Sample)){
  fq<-fastq_length[fastq_length$Sample==sample,]
  fs<-final[final$Sample==sample,]
  g<-ggplot()+ 
    geom_area(data=fq, aes(x=Length,y=RPM),fill="gray85",color="black")+
    geom_bar(data=fs, aes(x=Length,y=RPM, fill=Category),
             stat="identity", width=0.8, color=NA,size=0.73)+
    scale_fill_manual(values=allcolors)+
    theme_bw3() + 
    labs(x= "Read Length", 
         y="Total Reads Per Million",
         title=sample)+
    theme(plot.title = textTitle,
          axis.title = text20Bold,
          axis.text = text20,
          axis.line = element_line(colour = "gray75", size =0.73, linetype = "solid"),
          axis.ticks = element_line(size=0.73),axis.ticks.length=unit(0.3,"cm"),
          legend.text = text20Bold,
          legend.title=element_blank())
  print(g)  
}

dev.off()

cellWidth=1500
scales="free"

sampleCount<-length(unique(fastq_length$Sample))
facetColCount = ceiling(sqrt(sampleCount))
facetRowCount = ceiling( sampleCount * 1.0 / facetColCount)
width=max(2000, cellWidth * facetColCount + 300)
height=max(2000, cellWidth * facetRowCount)

png(file=paste0(outFile, ".png"), width=width, height=height, res=300)

g<-ggplot()+ 
  geom_area(data=fastq_length, aes(x=Length,y=RPM),fill="gray85",color="black")+
  geom_bar(data=final, aes(x=Length,y=RPM, fill=Category),
           stat="identity", width=0.8, color=NA,size=0.73)+
  scale_fill_manual(values=allcolors)+
  theme_bw3() + 
  labs(x = "Read Length", 
       y = "Total Reads Per Million")+
  theme(plot.title = textTitle,
        axis.title = text20Bold,
        axis.text = text20,
        axis.line = element_line(colour = "gray75", size=0.73, linetype = "solid"),
        axis.ticks = element_line(size=0.73),axis.ticks.length=unit(0.3,"cm"),
        strip.text = text20Bold,
        legend.text = text20Bold,
        legend.title = element_blank()) +
        facet_wrap(~Sample, ncol=facetColCount, scales=scales)
print(g)

dev.off()

groupFileList<-parSampleFile2
groups<-read.table(groupFileList, sep="\t", stringsAsFactor=F)
colnames(groups)<-c("Sample", "Group")
if(length(unique(groups$Group)) > 1){
  uniqueSamples<-unique(final$Sample)
  otherSamples<-uniqueSamples[!(uniqueSamples %in% groups$Sample)]
  if(length(otherSamples) > 0){
    groups<-rbind(groups, data.frame("Sample"=otherSamples, "Group"="Other"))
  }
  sampleMap<-split(groups$Group, groups$Sample)
  sample_group=lapply(sampleMap, function(x){paste(x, collapse="/")})

  final$Group=unlist(sample_group[as.character(final$Sample)])

  #saveRDS(fastq_length, file="fastq_length.rds")
  #saveRDS(final, file="final.rds")


  #fastq_length<-readRDS("c:/projects/temp/fastq_length.rds")
  sr<-split(fastq_length$RPM, fastq_length$Sample)
  sr_max<-lapply(sr, max)

  #final<-readRDS("c:/projects/temp/final.rds")

  fl<-unlist(sr_max[final$Sample])
  final$NRPM<-final$RPM / fl

  #allcolors<-c("blue","green","red", "brown","purple", "yellow", "black")[1:length(unique(final$Category))]

  #textTitle<-element_text(face= "bold", color = "black", size=22, hjust=0.5)
  #text20Bold<-element_text(face= "bold", color = "black", size=20)
  #text20<-element_text(color = "black", size=20)

  ## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
  ##   data: a data frame.
  ##   measurevar: the name of a column that contains the variable to be summariezed
  ##   groupvars: a vector containing names of columns that contain grouping variables
  ##   na.rm: a boolean that indicates whether to ignore NA's
  ##   conf.interval: the percent range of the confidence interval (default is 95%)
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE) {
      library(plyr)

      # New version of length which can handle NA's: if na.rm==T, don't count them
      length2 <- function (x, na.rm=FALSE) {
          if (na.rm) sum(!is.na(x))
          else       length(x)
      }

      # This does the summary. For each group's data frame, return a vector with
      # N, mean, and sd
      datac <- ddply(data, groupvars, .drop=.drop,
        .fun = function(xx, col) {
          c(N    = length2(xx[[col]], na.rm=na.rm),
            mean = mean   (xx[[col]], na.rm=na.rm),
            sd   = sd     (xx[[col]], na.rm=na.rm)
          )
        },
        measurevar
      )

      # Rename the "mean" column    
      datac <- rename(datac, c("mean" = measurevar))

      datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

      # Confidence interval multiplier for standard error
      # Calculate t-statistic for confidence interval: 
      # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
      ciMult <- qt(conf.interval/2 + .5, datac$N-1)
      datac$ci <- datac$se * ciMult

      return(datac)
  }

  slimse<-summarySE(final, measurevar="NRPM", groupvars=c("Length","Category", "Group"))
  #saveRDS(slimse, file="slimse.rds")
  n_group=length(unique(slimse$Group))
  png(paste0(outFile, ".group_category.png"), width=2000, height=400 * n_group + 100, res=300)
  gga<-ggplot(slimse, aes(x=Length, y=NRPM, color=Category)) + xlab("") + ylab("Relative RPM") + 
    geom_errorbar(aes(ymin=NRPM-se, ymax=NRPM+se), width=.5) +
    geom_line() +
    geom_point() +
    scale_color_manual(values=allcolors) +
    facet_grid(Group~.) + theme_bw() + xlim(15,40) + 
    theme_bw3()
  print(gga)
  dev.off()
}
