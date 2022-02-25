library(ggplot2)
library(data.table)
library(dplyr)

if(! exists("inputFile")){
  args <- commandArgs(TRUE)
  if(length(args) == 0){
    setwd("/scratch/cqs/shengq2/ravi_shah_projects/20220107_rnaseq_discovery_hg38/figures")
    inputFile<-"CD7.depth.txt"
    name<-'CD7'
    locus<-"chr17:82314868-82317602"
    width<-1000
    height<-800
  }else{
    inputFile<-args[1]
    name<-args[2]
    locus<-args[3]
    width<-as.numeric(args[4])
    height<-as.numeric(args[5])
  }
}

cat(inputFile, "\n")
dataForPlot<-fread(file=inputFile, sep="\t", header=T)
dataForPlot$FileName<-gsub("_", " ", dataForPlot$File)

outputFile = paste0(name, ".coverage.png")

png(outputFile, height=height, width=width, res=300)
g<-ggplot(dataForPlot, aes(Position, File)) + 
  geom_tile(aes(fill = PositionCount))+
  scale_fill_gradientn(colours=rev(heat.colors(100)))+   
  ylab("") + xlab(paste0(name, " [", locus, "]"))+
  theme_bw() +  
  theme(strip.background=element_blank())
print(g)
dev.off()
