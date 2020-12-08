library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  inputFile = "/scratch/cqs/alexander_gelbard_projects/20201202_5126_scRNA_split/split_samples_summary/result/scRNA_5126__fileList1.list"
  outputPrefix = "/scratch/cqs/alexander_gelbard_projects/20201202_5126_scRNA_split/split_samples_summary/result/scRNA_5126.HTO.summary"
}else{
  inputFile = args[1]
  outputPrefix = args[2]
}

files=read.table(inputFile, sep="\t", stringsAsFactor=F)

dat=apply(files, 1, function(x){
  dname=x[[2]]
  dfile=x[[1]]
  dd=read.csv(dfile)
  table(dd$HTO.global)
})
colnames(dat)=files$V2
mdat=melt(dat)
colnames(mdat)=c("Class", "Sample", "Cell")

pdf(paste0(outputPrefix, ".global.pdf"))
g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw()
print(g)
dev.off()

dfile=files$V1[1]
dat=apply(files, 1, function(x){
  dname=x[[2]]
  dfile=x[[1]]
  dd=read.csv(dfile, stringsAsFactor=F)
  dd$Class=dd$HTO.global
  dd$Class[dd$Class=="Singlet"] = dd$HTO[dd$Class=="Singlet"]
  df=data.frame(table(dd$Class))
  df$Sample=dname
  df
})
mdat=do.call(rbind, dat)
colnames(mdat)=c("Class", "Cell", "Sample")
mdat$Class=factor(mdat$Class, levels=sort(unique(as.character(mdat$Class))))

pdf(paste0(outputPrefix, ".hto.pdf"))
g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw()
print(g)
dev.off()
