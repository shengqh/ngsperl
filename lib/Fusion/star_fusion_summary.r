rm(list=ls()) 
outFile='rnaseq'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/jennifer_pietenpol_projects/20241213_CALGB40603_rnaseq_hg38/star_fusion_summary/result')

### Parameter setting end ###

library(data.table)
library(dplyr)
library(ggplot2)

filemap=read.table(parSampleFile1, header=F)

fusion_count=unlist(lapply(filemap$V1, function(x){
  filedata=fread(x)
  return(nrow(filedata))
}))

fusion_df=data.frame("Sample"=filemap$V2, "FusionCount"=fusion_count)
write.csv(fusion_df, file=paste0(outFile, ".fusion_count.csv"), row.names=F)

g=ggplot(fusion_df, aes(x="Fusion", y=FusionCount)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0(outFile, ".fusion_count.png"), g, width=2, height=2, units="in", dpi=300, bg="white")

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')
