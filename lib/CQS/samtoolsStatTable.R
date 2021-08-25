options(bitmapType='cairo')

resultPrefix<-outFile
statFileList<-parSampleFile1

library(ggplot2)
library(reshape2)
library(stringr)

#resultPrefix<-"/workspace/shengq1/guoyan/20160922_rnaediting/fastqc_count_vis/result/stat"
#statFileList<-"/workspace/shengq1/guoyan/20160922_rnaediting/fastqc_count_vis/result/fileList1.txt"

statFiles=read.delim(statFileList, header=F, stringsAsFactor=F)

tbl=apply(statFiles, 1, function(x){
  file=x[1]
  statName=x[2]
  
  stat=read.delim(file, header=F,check.names=F)
  stat=head(stat, -1)
  
  reads=data.frame(Reads=str_match(stat$V1, "\\d+"))
  colnames(reads)=statName
  names=str_match(stat$V1, "\\d+\\s*\\+\\s*\\d+\\s*([^(]*)")[,2]
  names=gsub("\\s+$", "", names)
  rownames(reads)=names
  reads
})

tbl2=do.call("cbind", tbl)
tbl3=data.frame(t(tbl2))
tbl4=cbind(data.frame(Sample=rownames(tbl3)), tbl3)
write.table(tbl4,paste0(resultPrefix,".Reads.tsv"),quote=F,row.names=F,sep="\t")

