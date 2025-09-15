rm(list=ls()) 
outFile='RA_4893_2'
parSampleFile1=''
parSampleFile2=''
parSampleFile3=''
parFile1='/scratch/vickers_lab/projects/20220707_4893_2_RA_smRNA_mouse_v5_forSpcount/nonhost_genome/bowtie1_bacteria_group1_pm_table/result/bacteria_group1_pm_RA_4893_2.read.count'
parFile2='/scratch/vickers_lab/projects/20220707_4893_2_RA_smRNA_mouse_v5_forSpcount/nonhost_genome/bowtie1_bacteria_group2_pm_table/result/bacteria_group2_pm_RA_4893_2.read.count'
parFile3='/scratch/vickers_lab/projects/20220707_4893_2_RA_smRNA_mouse_v5_forSpcount/nonhost_genome/refseq_bacteria_table/result/RA_4893_2.read.count'


setwd('/scratch/vickers_lab/projects/20220707_4893_2_RA_smRNA_mouse_v5_forSpcount/data_visualization/bacteria_count_summary/result')

### Parameter setting end ###

source("countTableVisFunctions.R")

options(bitmapType='cairo')

library(ggplot2)
library(reshape2)

group1<-read.table(parFile1, sep="\t", stringsAsFactors = F, header=T, row.names=1)
group2<-read.table(parFile2, sep="\t", stringsAsFactors = F, header=T, row.names=1)
refseq<-read.table(parFile3, sep="\t", stringsAsFactors = F, header=T, row.names=1)

gs1<-colSums(group1)
gs2<-colSums(group2)
rs<-colSums(refseq)

g12<-rbind(group1, group2[!(rownames(group2) %in% rownames(group1)),,drop=F])
gs12<-colSums(g12)

df<-data.frame("Microbiome"=gs1, "Environment"=gs2, "Microbiome_Enviroment"=gs12, "Refseq"=rs)
df<-t(df)

write.csv(df, paste0(outFile, ".csv"))

mdf<-melt(df)
colnames(mdf)<-c("Category", "Sample", "Reads")

g<-ggplot(mdf, aes(x=Category, y=Reads)) + 
  geom_bar(stat="identity") + 
  facet_grid(rows=Sample~.) + 
  coord_flip() + 
  theme_bw3()

width=2000
height=min(10000, max(2000, ncol(df) * 400))

png(file=paste0(outFile, ".bar1.png"), width=width, height=height, res=300)
print(g)
dev.off()

g2_a<-group2[!(rownames(group2) %in% rownames(group1)),]
refseq_a<-refseq[!(rownames(refseq) %in% c(rownames(group1), rownames(g2_a))),]

gs2_a<-colSums(g2_a)
rs_a<-colSums(refseq_a)

df<-data.frame("Microbiome"=gs1, "Environment+"=gs2_a, "Refseq+"=rs_a, check.names=F)
df<-t(df)

mdf<-melt(df)
colnames(mdf)<-c("Category", "Sample", "Reads")
mdf$Category<-factor(mdf$Category, levels=c("Refseq+", "Environment+", "Microbiome"))

write.csv(mdf, paste0(outFile, ".2.csv"))

width=3000
height=2000

g<-ggplot(mdf, aes(x=Sample, y=Reads, fill=Category)) + geom_bar(stat="identity") + theme_bw3() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
png(file=paste0(outFile, ".bar2.png"), width=width, height=height, res=300)
print(g)
dev.off()

temp<-tapply(mdf[,"Reads"],mdf[,"Sample"],sum)
mdf$Percentage<-mdf[,"Reads"]/temp[mdf[,"Sample"]] * 100

write.csv(mdf, paste0(outFile, ".3.csv"))
g<-ggplot(mdf, aes(x=Sample, y=Percentage, fill=Category)) + geom_bar(stat="identity") + theme_bw3() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
png(file=paste0(outFile, ".bar3.png"), width=width, height=height, res=300)
print(g)
dev.off()
