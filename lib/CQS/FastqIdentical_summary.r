rm(list=ls()) 
outFile='apd_smallrna_hg38'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/cqs/ravi_shah_projects/shengq2/20230201_apd_smallrna_hg38/preprocessing/fastqc_post_trim_summary/result/apd_smallrna_hg38.FastQC.reads.tsv'
parFile2=''
parFile3=''


setwd('/scratch/cqs/ravi_shah_projects/shengq2/20230201_apd_smallrna_hg38/preprocessing/identical_summary/result')

### Parameter setting end ###
library(data.table)
library(ggplot2)

files_tb<-read.table(parSampleFile1, sep="\t", header=F)
files_map<-unlist(split(files_tb$V1, files_tb$V2))

total_reads<-read.table(parFile1, sep="\t", header=T)
total_reads_map<-unlist(split(total_reads$Reads, total_reads$Sample))

all_df = NULL
sample = names(files_map)[1]
for(sample in names(files_map)){
  dup_file = files_map[sample]
  cat(sample, ":", dup_file, "\n")
  dup_count = fread(dup_file, nrows=10000)
  total_read = total_reads_map[sample]
  dup_count$log_rpm = log(dup_count$Count / total_read * 1000000)
  dup_count$seq_length = nchar(dup_count$Sequence)
  dup_count$sample = sample
  dup_count = dup_count[,c("sample", "log_rpm", "seq_length")]
  all_df = rbind(all_df, dup_count)
}

#library(MASS)
#bandwidth.nrd(all_df$seq_length)
#bandwidth.nrd(all_df$log_rpm)

nsample=length(names(files_map))
nwidth=ceiling(sqrt(nsample))
nheight=ceiling(nsample / nwidth)

png("density_2d.png", width=500 * nwidth + 200, height=500 * nheight + 100, res=300)
g<-ggplot(all_df, aes(x=seq_length, y=log_rpm)) + geom_density_2d_filled(contour_var = "ndensity") + facet_wrap(~sample) + theme_bw()
print(g)
dev.off()

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
