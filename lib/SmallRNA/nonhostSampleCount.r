rm(list=ls()) 
outFile='CM_8643_bakeoff'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/scratch/vickers_lab/projects/20220830_8643_CM_smRNA_human_bakeoff/nonhost_genome/nonhost_genome_count/result')

### Parameter setting end ###

library(data.table)

files<-read.table(parSampleFile1, header=F, sep="\t")

myoptions<-read.table(parSampleFile2, header=F, sep="\t")
myoptions<-split(myoptions$V1, myoptions$V2)

all_reads<-NULL
file<-files$V1[1]
for(file in files$V1){
  cat("reading ", file, "...\n")
  cur_reads<-fread(file, sep="\t", header=T)
  all_reads<-rbind(all_reads, cur_reads)
}

unique_reads<-all_reads[!duplicated(all_reads$Sequence), c(2:ncol(all_reads)), with=FALSE]
rm(all_reads)

sample_reads<-apply(unique_reads, 2, sum)

sf<-data.frame("Sample" = names(sample_reads), "Count"=sample_reads)
write.table(sf, paste0(outFile, myoptions$extension), sep="\t", quote=F, col.names=T, row.names=F)
