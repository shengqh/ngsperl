
library(ggplot2)
library(reshape2)

inputFile=parSampleFile1
nameMapFile=parSampleFile2
outputPrefix=paste0(outFile, ".HTO.summary")

cat("inputFile=", inputFile, "\n")
cat("outputPrefix=", outputPrefix, "\n")
cat("nameMapFile=", nameMapFile, "\n")

files=read.table(inputFile, sep="\t", stringsAsFactor=F)

dat=apply(files, 1, function(x){
  dname=x[[2]]
  dfile=x[[1]]
  dd=read.csv(dfile, check.names=F)
  table(dd$HTO.global)
})
colnames(dat)=files$V2
write.csv(dat, file=paste0(outputPrefix, ".csv"), quote=F)

mdat=reshape2::melt(dat)
colnames(mdat)=c("Class", "Sample", "Cell")

png(paste0(outputPrefix, ".global.png"), width=1600, height=1200, res=300)
g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust=0.5)) + xlab("")
print(g)
dev.off()

hasNameMap = (!is.na(nameMapFile)) && (nameMapFile != "")
if (hasNameMap) {
  df = read.table(nameMapFile, sep="\t", stringsAsFactor=F)
  namemap=split(df$V1, df$V2)
}

dat=apply(files, 1, function(x){
  dname=x[[2]]
  dfile=x[[1]]
  dd=read.csv(dfile, stringsAsFactors=F, check.names=F)
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
  mdat$Cell[mdat$Var1 %in% names(namemap)] = unlist(namemap[mdat$Var1[mdat$Var1 %in% names(namemap)]])
}
mdat=mdat[,c("Sample", "Cell", "Freq")]
colnames(mdat)=c("Sample", "Class", "Cell")

dat=reshape2::dcast(mdat, formula="Sample~Class")
write.csv(dat, file=paste0(outputPrefix, ".csv"), quote=F, row.names=F)

png(paste0(outputPrefix, ".sample.png"), width=1600, height=1200, res=300)
g<-ggplot(mdat, aes(x=Sample, y=Cell, fill=Class, label=Cell)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 0, vjust=0.5)) + xlab("")
print(g)
dev.off()
