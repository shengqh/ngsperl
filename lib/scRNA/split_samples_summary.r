library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  inputFile = "/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/hto_samples_cutoff_summary/result/scRNA__fileList1.list"
  nameMapFile = "/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/hto_samples_cutoff_summary/result/scRNA__fileList2.list"
  outputPrefix = "/scratch/cqs/shengq2/paula_hurley_projects/20210303_scRNA_human/hto_samples_cutoff_summary/result/scRNA.HTO.summary"
}else{
  inputFile = args[1]
  outputPrefix = args[2]
  nameMapFile = args[3]
}
cat("inputFile=", inputFile, "\n")
cat("outputPrefix=", outputPrefix, "\n")
cat("nameMapFile=", nameMapFile, "\n")

files=read.table(inputFile, sep="\t", stringsAsFactor=F)

dat=apply(files, 1, function(x){
  dname=x[[2]]
  dfile=x[[1]]
  dd=read.csv(dfile)
  table(dd$HTO.global)
})
colnames(dat)=files$V2
write.csv(dat, file=paste0(outputPrefix, ".csv"), quote=F)

mdat=melt(dat)
colnames(mdat)=c("Class", "Sample", "Cell")

pdf(paste0(outputPrefix, ".global.pdf"))
g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw()
print(g)
dev.off()

hasNameMap = ! is.na(nameMapFile)
if (hasNameMap) {
  df = read.table(nameMapFile, sep="\t", stringsAsFactor=F)
  namemap=df$V1
  names(namemap)=df$V2
}

dat=apply(files, 1, function(x){
  dname=x[[2]]
  dfile=x[[1]]
  dd=read.csv(dfile, stringsAsFactors=F)
  dd$HTO.final=dd$HTO.global
  dd$HTO.final[dd$HTO.global=="Singlet"] = dd$HTO[dd$HTO.global=="Singlet"]
  dt=data.frame(table(dd$HTO.final), stringsAsFactors=F)
  dt$Var1=as.character(dt$Var1)
  dt$Sample=dname
  dt
})
mdat=do.call("rbind", dat)
mdat$Cell=mdat$Var1
if (hasNameMap) {
  mdat$Cell[mdat$Var1 %in% names(namemap)] = namemap[mdat$Var1[mdat$Var1 %in% names(namemap)]]
}
mdat=mdat[,c("Sample", "Cell", "Freq")]
colnames(mdat)=c("Sample", "Class", "Cell")

dat=dcast(mdat, formula="Sample~Class")
write.csv(dat, file=paste0(outputPrefix, ".csv"), quote=F, row.names=F)

pdf(paste0(outputPrefix, ".global.pdf"))
g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + geom_bar(position="stack", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw()
print(g)
dev.off()
