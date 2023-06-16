
library(ggplot2)
library(reshape2)

microbial<-read.delim(parFile2, header=T, row.names=1,check.names=F)

all<-read.csv(parFile1, row.names=1,check.names=F)
all<-data.frame(t(all))
all$all<-rowSums(all)
all$Host<-(all$Host.Small.RNA + all$Mapped.to.Host.Genome) / all$all * 100
all$Microbial<-(microbial[rownames(all), "Count"]) / all$all * 100
all$Other<-100 - all$Host - all$Microbial

figData<-all[,c("Host", "Microbial", "Other")]

hasGroup = parSampleFile1 != ""
if(hasGroup){
  group<-read.delim(parSampleFile1, stringsAsFactors = F, header=F)
  if (exists("uniqueGroupNames")){
    group<-group[group$V2 %in% uniqueGroupNames,]
  }
  rownames(group)<-group$V1
  figData$Group<-group[rownames(figData), "V2"]
} else {
  figData$Group<-"All"
}
figData$Sample<-rownames(figData)

figData<-figData[!(rownames(figData) %in% c("B020", "B114", "D031")),]
figData<-figData[order(figData$Microbial, decreasing = T),]

mFigData<-melt(figData, id.vars=c("Sample", "Group"))
colnames(mFigData)<-c("Sample", "Group", "Category", "Percentage")
mFigData$Sample<-factor(mFigData$Sample, levels=figData$Sample)
mFigData$Category<-factor(mFigData$Category, levels=c("Microbial", "Host", "Other"))

mFigData$Group[is.na(mFigData$Group)] = "Unknown group"
write.table(mFigData, file=paste0(outFile, ".txt"), sep="\t", row.names=F, col.names=T, quote=F)

pdf(paste0(outFile, ".pdf"), width=10, height=6)
colors<-c("Microbial" = "chartreuse3", "Host" = "deepskyblue", "Other" = "gray")
g<-ggplot(mFigData) + geom_bar(aes(y = Percentage, x = Sample, fill = Category), stat="identity", width=1) + facet_grid(~Group, scales = "free_x") +
  scale_fill_manual(values=colors) +
  ylab("Percentage of Reads") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background =element_blank(), 
        panel.background = element_blank(),
        strip.text.x.top = element_text(angle = 90, hjust = 0))
print(g)
dev.off()
