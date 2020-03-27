setwd("/data/stein_lab/mjo_sRNA_data/20190717_michelle_smallRNA_2868_human_paper_figures")
library(ggplot2)
library(reshape2)

microbial<-read.delim("2868.microbial.tsv", header=T, row.names=1)

all<-read.csv("/data/stein_lab/mjo_sRNA_data/20170206_michelle_smallRNA_2868_human/data_visualization/reads_in_tasks/result/2868.NonParallel.TaskReads.csv",row.names=1)
all<-data.frame(t(all))
all$all<-rowSums(all)
all$Human<-(all$Host.Small.RNA + all$Mapped.to.Host.Genome) / all$all * 100
all$Microbial<-(microbial[rownames(all), "Count"]) / all$all * 100
all$Other<-100 - all$Human - all$Microbial

figData<-all[,c("Human", "Microbial", "Other")]

group<-read.delim("/data/stein_lab/mjo_sRNA_data/20170206_michelle_smallRNA_2868_human/data_visualization/reads_in_tasks/result/fileList2_pie.txt", stringsAsFactors = F, header=F)
group<-group[group$V2 %in% c("RA", "Control"),]
rownames(group)<-group$V1
figData$Group<-group[rownames(figData), "V2"]
figData$Sample<-rownames(figData)

figData<-figData[!(rownames(figData) %in% c("B020", "B114", "D031")),]
figData<-figData[order(figData$Microbial, decreasing = T),]

mFigData<-melt(figData, id.vars=c("Sample", "Group"))
colnames(mFigData)<-c("Sample", "Group", "Category", "Percentage")
mFigData$Sample<-factor(mFigData$Sample, levels=figData$Sample)
mFigData$Category<-factor(mFigData$Category, levels=c("Microbial", "Human", "Other"))
mFigData$Group<-factor(mFigData$Group, levels=c("RA", "Control"))

pdf("2868_human_reads_bar.pdf", width=10, height=6)
colors<-c("Microbial" = "chartreuse3", "Human" = "deepskyblue", "Other" = "gray")
g<-ggplot(mFigData) + geom_bar(aes(y = Percentage, x = Sample, fill = Category), stat="identity") + facet_grid(~Group, scales = "free_x") +
  scale_fill_manual(values=colors) +
  ylab("Percentage of Reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(g)
dev.off()
