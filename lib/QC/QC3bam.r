# TODO: Add comment
#
# Author: zhaos
###############################################################################

args = commandArgs(trailingOnly=TRUE)

resultDir<-args[1]
outFile<-args[2]

setwd(resultDir)

library(reshape2)
library(ggplot2)
library(patchwork)
library(dplyr)

fileName<-paste(resultDir,'/bamResult/bamSummary.txt',sep="")
if (!file.exists(fileName)) {
	fileName<-paste(getwd(),'/',fileName,sep="")
}
allResult<-read.delim(fileName,header=T,row.names=1,check.names=F)

figureDir<-paste(resultDir,"/bamFigure/",sep="")
dir.create(figureDir, showWarnings = FALSE)
oldwd<-getwd()
setwd(figureDir)
file.remove(list.files("."))
if (ncol(allResult)>=35+1) {
	colList<-lapply(list(5:11,c(12:14,NA,15:17),c(18:20,NA,21:23),c(24:26,NA,27:29),c(30:32,NA,33:35)),function(x){x+1})
} else if (ncol(allResult)>=29+1) {
	colList<-lapply(list(5:11,c(12:14,NA,15:17),c(18:20,NA,21:23),c(24:26,NA,27:29)),function(x){x+1})
} else {
	colList<-lapply(list(5:11,c(12:14,NA,15:17),c(18:20,NA,21:23)),function(x){x+1})
}
if (ncol(allResult)>=35+1) {
	colBox<-(5:35)+1
} else if (ncol(allResult)>=29+1) {
	colBox<-(5:29)+1
} else {
	colBox<-(5:ncol(allResult))+1
}

slimResult=allResult[,c(1,6,8,10:12)]
if(all(slimResult$Unmapped==0)) {
  slimResult$Unmapped<-NULL
}
if(all(slimResult$`Off-target-mito`==0)) {
  slimResult$`Off-target-mito`<-NULL
}

mres=reshape2::melt(slimResult, id.var="SM")
colnames(mres)=c("Sample", "Category", "Reads")

mres$Category=dplyr::recode(mres$Category, `On-target`="Exon", `Off-target-intron`="Intron", `Off-target-intergenic`="Intergenic")

g1=ggplot(mres, aes(x=Sample, y=Reads, fill=Category)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title="", x="", y="Reads") + scale_fill_brewer(palette="Set1")

g2=ggplot(mres, aes(x=Sample, y=Reads, fill=Category)) + 
  geom_bar(position="fill", stat="identity") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title="", x="", y="Percent") + scale_fill_brewer(palette="Set1")

g=g1/g2
ggsave(paste0(outFile, ".readsByCategory.png"), g, width=6, height=6, units="in")
