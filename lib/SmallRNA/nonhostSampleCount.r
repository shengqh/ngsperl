rm(list=ls()) 
outFile='GSE137617'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/scratch/cqs/shengq2/bugfix/nonhost_count')

### Parameter setting end ###

library(data.table)

files<-read.table(parSampleFile1, header=F, sep="\t")

myoptions<-read.table(parSampleFile2, header=F, sep="\t")
myoptions<-split(myoptions$V1, myoptions$V2)

all_reads<-NULL
file<-files$V1[6]
for(file in files$V1){
  cat("reading ", file, "...\n")
  cur_reads<-fread(file, sep="\t", header=T)
  if (is.null(all_reads)){
    all_reads = cur_reads
  }else{
    if(ncol(all_reads) != ncol(cur_reads)){
      all_colnames = sort(unique(c(colnames(all_reads)[2:ncol(all_reads)], colnames(cur_reads)[2:ncol(cur_reads)])))
      for(colname in all_colnames){
        if(!(colname %in% colnames(all_reads))){
          all_reads[,colname] = 0
        }
        if(!(colname %in% colnames(cur_reads))){
          cur_reads[,colname] = 0
        }
      }
      all_reads<-all_reads[,c("Sequence", all_colnames),with=FALSE]
      cur_reads<-cur_reads[,c("Sequence", all_colnames),with=FALSE]
    }
    all_reads<-rbind(all_reads, cur_reads)
  }
}

unique_reads<-all_reads[!duplicated(all_reads$Sequence), c(2:ncol(all_reads)), with=FALSE]
rm(all_reads)

sample_reads<-apply(unique_reads, 2, sum)

sf<-data.frame("Sample" = names(sample_reads), "Count"=sample_reads)
write.table(sf, paste0(outFile, myoptions$extension), sep="\t", quote=F, col.names=T, row.names=F)
