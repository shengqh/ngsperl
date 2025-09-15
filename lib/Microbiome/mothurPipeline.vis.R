
#Instructions####
#https://rpubs.com/dillmcfarlan/R_microbiotaSOP


sampleGroupFile=parSampleFile1
mothur_shared_file=parFile1
mothur_constaxonomy_file=parFile2
totalReads_file=parFile3

#dataDir="/scratch/cqs/zhaos/WilsonKeith/20210309_5932AG16S/pipeline/mothur_pipeline/result"
#mothur_shared_file="stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
#mothur_constaxonomy_file="stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"

#workDir="/scratch/cqs/zhaos/WilsonKeith/20210309_5932AG16S/pipeline/mothur_pipeline_vis/result"


###################################1
#prepare data by phyloseq####
###################################1
library(phyloseq)
#setwd(dataDir)
pseqAll=import_mothur(mothur_shared_file=mothur_shared_file,mothur_constaxonomy_file=mothur_constaxonomy_file)

#setwd(workDir)

#remove Mock or control samples
#expSamples=setdiff(sample_names(pseqAll),c("AG.NC","AG.PC"))
#pseq <- prune_samples(expSamples, pseqAll)
pseq=pseqAll

#sample information
sampleGroups=read.delim(sampleGroupFile,header=FALSE,as.is=FALSE)
colnames(sampleGroups)=c("Sample","group")
row.names(sampleGroups)=sampleGroups$Sample
SDF = sampleGroups[sample_names(pseq),-1,drop=FALSE]
#SDF
sample_data(pseq) <- sample_data(SDF)


###################################1
#total reads and aligned reads####
###################################1

if (totalReads_file !="" & file.exists(totalReads_file)) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  readsOutFile=paste0(outFile,".reads")
  
  totalReadsTable=read.delim(totalReads_file,header=TRUE,as.is=TRUE)
  totalReadsTable=unique(totalReadsTable[,c("Sample","Reads")])
  totalReads=totalReadsTable$Reads
  names(totalReads)=totalReadsTable$Sample
  
  otuTotalReads=colSums(otu_table(pseq))
  
  dataForPlot=data.frame(Sample=row.names(sample_data(pseq)),sample_data(pseq),
                         #TotalReads=totalReads[row.names(sample_data(pseq))],
                         ClassfiedReads=otuTotalReads[row.names(sample_data(pseq))])
  dataForPlot$NotClassfiedReads=totalReads[dataForPlot$Sample]-dataForPlot$ClassfiedReads
  
  dataForPlot=dataForPlot %>% pivot_longer(cols=c("ClassfiedReads","NotClassfiedReads"))
  
  p=ggplot(dataForPlot,aes(x=Sample,y=value,fill=name))+
    facet_grid(rows=~group,scales = "free_x",space="free")+ylab("Reads")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p1=p+geom_bar(stat = "identity")
  p2=p+geom_bar(stat = "identity",position="fill")+ scale_y_continuous(labels = scales::percent_format())

  pdf(paste0(readsOutFile,".pdf"),width=9)
  plot(p1)
  dev.off()
  
  pdf(paste0(readsOutFile,".percent.pdf"),width=9)
  plot(p2)
  dev.off()
}



###################################1
#use microbiome for analysis####
###################################1

library(microbiome)

## Alpha diversity####
tab1=estimate_richness(pseq)
groupMean=apply(tab1[,c("Chao1","ACE","Shannon","Simpson")],2,function(x) tapply(x,meta(pseq)$group,mean))
groupSD=apply(tab1[,c("Chao1","ACE","Shannon","Simpson")],2,function(x) tapply(x,meta(pseq)$group,sd))

resultOut=paste(round(groupMean,2),"+",round(groupSD,2))
resultOut=matrix(resultOut,ncol=4)
colnames(resultOut)=colnames(groupMean)
row.names(resultOut)=row.names(groupMean)

AlphaDiversityOutFile=paste0(outFile,".AlphaDiversity")

write.csv(tab1,paste0(AlphaDiversityOutFile,".csv"))
write.csv(resultOut,paste0(AlphaDiversityOutFile,".GroupSummary.csv"))

###################################1
## Microbiome composition####
###################################1
library(tidyverse)
library(ggplot2)

pseqCompositionalAll <- transform(pseq, "compositional")
#rank_names(pseqCompositionalAll)
#"Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"

for (level in c("Rank2","Rank5","Rank6")) {
#  level="Rank5"
  pseqCompositional <- aggregate_rare(pseqCompositionalAll, level = level, detection = 1/100, prevalence = 10/100)
  
  compositionOutFile=paste0(outFile,".compositionOfTax",level)
  
  #Composition barplots
  p <- pseqCompositional %>%
    plot_composition(otu.sort = "abundance",group_by="group") +
    scale_y_continuous(label = scales::percent)
  pdf(paste0(compositionOutFile,".pdf"),width=9)
  print(p)
  dev.off()
  write.csv(otu_table(pseqCompositional),paste0(compositionOutFile,".csv"))
  
  #average by group
  p <- pseqCompositional %>% plot_composition(otu.sort = "abundance", average_by = "group")+ scale_y_continuous(label = scales::percent)
  pdf(paste0(compositionOutFile,".groupAverage.pdf"),width=9)
  print(p)
  dev.off()
  #export group level data
  pseqCompositionalOtuGroupMean=apply(otu_table(pseqCompositional),1,function(x) tapply(x,meta(pseqCompositional)$group,mean))
  pseqCompositionalOtuGroupMean=t(pseqCompositionalOtuGroupMean)
  write.csv(pseqCompositionalOtuGroupMean,paste0(compositionOutFile,".groupAverage.csv"))
  
  #kruskal.test
  percentData=otu_table(pseqCompositional)
  groupData=meta(pseq)$group
  resultPAll=NULL
  for (i in 1:nrow(percentData)) {
    #i=1
    resultP=kruskal.test(as.vector(percentData[i,])~ groupData)$p.value
    resultPAll=c(resultPAll,resultP)
  }
  resultOut=cbind(ID=row.names(percentData),p=resultPAll,adjustedP=p.adjust(resultPAll,method="fdr"))
  write.csv(resultOut,paste0(compositionOutFile,".testResult.csv"),row.names=FALSE)
  
  
  #Composition heatmaps
  p <- plot_composition(pseqCompositional,plot.type = "heatmap",sample.sort = "neatmap", otu.sort = "neatmap")
  pdf(paste0(compositionOutFile,".heatmap.pdf"),width=9)
  print(p+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  dev.off()

  
  # #Plot taxa prevalence
  # p0 <- core(pseq, detection = 0.1/100, prevalence = 1/100)
  # plot_taxa_prevalence(p0, level, detection = 0.1/100)
}

# level="Rank2"
# pseqCompositional <- transform(pseq, "compositional")
# pseqCompositional <- aggregate_rare(pseqCompositional, level = level, detection = 0.1/100, prevalence = 10/100)


###################################1
## PERMANOVA####
###################################1
library(vegan)
PERMANOVAOutFile=paste0(outFile,".PERMANOVA")

pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

#save output
sink(paste0(PERMANOVAOutFile,".model.txt"))

permanova <- adonis(t(otu) ~ group, data = meta, permutations=1000, method = "bray")
#permanova <- adonis(t(otu) ~ WtOrSmox*Spd, data = meta, permutations=1000, method = "bray")

# P-value
permanova
#print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])


# Investigate the top factors in PERMANOVA
coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
pdf(paste0(PERMANOVAOutFile,".TopCoef.pdf"))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top OTU")
dev.off()

emptyInd=which(apply(otu,1,function(x) all(x==0)))
if (length(emptyInd)>0) {
  otuNoEmpty=otu[-emptyInd,]
} else {
  otuNoEmpty=otu
}
BC.dist=vegdist(t(otuNoEmpty), distance="bray")
adonis(BC.dist ~ group, data = meta, permutations = 1000)
#adonis(BC.dist ~ WtOrSmox*Spd, data = meta, permutations = 1000)

#stop save output
sink()

#p <- plot_landscape(pseq.rel,col = "group", size = 3,method = "PCoA")
#p

gp_bray_pcoa = ordinate(pseq.rel, "PCoA", "bray")
#bugfix for only one annotation column in sample_data
if (ncol(sample_data(pseq.rel))==1) { #https://github.com/joey711/phyloseq/issues/541
  sample_data(pseq.rel)$Group1=sample_data(pseq.rel)$group
}
pdf(paste0(PERMANOVAOutFile,".PCA.pdf"))
p <- plot_ordination(pseq.rel,gp_bray_pcoa,col = "group")
plot(p+ stat_ellipse() +  theme_bw() + xlab("PCoA1")+ ylab("PCoA2"))
dev.off()

#export corrdinate to Alain so that he can make figure by himself
write.csv(gp_bray_pcoa$vectors[,1:2],paste0(PERMANOVAOutFile,".PCA.csv"))


save.image(paste0(outFile,".RData"))
