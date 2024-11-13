#renv::restore()
print("BEGIN CombineDECoN.R")

library(R.utils)
library(optparse)
library(data.table)
library(dplyr)

########## Process inputs #######

option_list<-list(
    make_option("--rdatas",help="text file containing list of DECoN count RData files",dest='rdatas'),
    make_option("--out",default="DECoN",help="output prefix, default =DECoN",dest='out')
)

opt<-parse_args(OptionParser(option_list=option_list))

rdata_file=opt$rdatas                                                                        #location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
output=opt$out                                                               #location and name of file to save the output to; will be saved as 'output_counts.RData'

DEBUG=FALSE
if(DEBUG){
  setwd('/nobackup/h_cqs/shengq2/temp/panel_exomeseq/T03_DECoN_makeCNVcalls/result')
  rdata_file="/nobackup/h_cqs/shengq2/temp/panel_exomeseq/T03_DECoN_makeCNVcalls/result/fileList1.list"
  output="panal_test"
}else{
  if(rdata_file=="NULL"){
    rdata_file=NULL
  }

  if(is.null(rdata_file)){
    print("ERROR rdata files must be provided. Execution halted")
    quit()
  }

  print(output)
}
#################################

print(rdata_file)

rdata_files=fread(rdata_file, header=FALSE)
if( nrow(rdata_files)==0){
  print("ERROR NO RDATA FILES DETECTED")
}

rdata_map=split(rdata_files$V1, rdata_files$V2)

all_sample_names<-names(rdata_map) 
names(all_sample_names)<-NULL

final_counts<-NULL
cur_sample = all_sample_names[1]
for(cur_sample in all_sample_names){
  rdata=rdata_map[[cur_sample]]
  print(paste("Reading ", rdata, " ..."))
  load(rdata)
  colnames(counts)[ncol(counts)]<-cur_sample
  if(is.null(final_counts)){
    final_counts<-counts
  }else{
    final_counts<-cbind(final_counts, counts[,ncol(counts),drop=FALSE])
  }
}

counts=final_counts
sample.names=all_sample_names
save(counts,bams,bed.file,sample.names,fasta,file=paste(output,".RData",sep=""))                                   #saves workspace as 'output_counts.RData'
