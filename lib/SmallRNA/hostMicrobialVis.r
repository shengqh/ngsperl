
library(ggplot2)
library(reshape2)

microbial<-read.delim(parFile2, header=T, row.names=1)

all<-read.csv(parFile1, row.names=1)
all<-data.frame(t(all))
all$all<-rowSums(all)
all$Host<-(all$Host.Small.RNA + all$Mapped.to.Host.Genome) / all$all * 100
all$Microbial<-(microbial[rownames(all), "Count"]) / all$all * 100
all$Other<-100 - all$Host - all$Microbial

figData<-all[,c("Host", "Microbial", "Other")]

group<-read.delim(parSampleFile1, stringsAsFactors = F, header=F)
if (exists("groupNames")){
  group<-group[group$V2 %in% groupNames,]
}
rownames(group)<-group$V1
figData$Group<-group[rownames(figData), "V2"]
figData$Sample<-rownames(figData)

figData<-figData[!(rownames(figData) %in% c("B020", "B114", "D031")),]
figData<-figData[order(figData$Microbial, decreasing = T),]

mFigData<-melt(figData, id.vars=c("Sample", "Group"))
colnames(mFigData)<-c("Sample", "Group", "Category", "Percentage")
mFigData$Sample<-factor(mFigData$Sample, levels=figData$Sample)
mFigData$Category<-factor(mFigData$Category, levels=c("Microbial", "Host", "Other"))

pdf(paste0(outFile, ".pdf"), width=10, height=6)
colors<-c("Microbial" = "chartreuse3", "Host" = "deepskyblue", "Other" = "gray")
g<-ggplot(mFigData) + geom_bar(aes(y = Percentage, x = Sample, fill = Category), stat="identity") + facet_grid(~Group, scales = "free_x") +
  scale_fill_manual(values=colors) +
  ylab("Percentage of Reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background =element_blank(), 
        panel.background = element_blank())
print(g)
dev.off()