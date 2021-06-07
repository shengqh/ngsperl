# rm(list=ls()) 
# outFile='output'
# parSampleFile1='fileList2.txt'
# parSampleFile2=""
# parSampleFile3=''
# parSampleFile4=''
# parFile1='RA_4949_mouse.Category.Table.csv'
# parFile2=''
# parFile3=''

#setwd("C:/projects/composition_test")

library(reshape2)
library(ggplot2)
library(DirichletReg)
library(pheatmap)

comp<-read.csv(parFile1,row.names=1, check.names=F)

getSampleInGroup<-function(groupDefineFile, samples, useLeastGroups=FALSE,onlySamplesInGroup=FALSE){
  allGroupData<-read.delim(groupDefineFile,as.is=T,header=F)
  if(ncol(allGroupData) < 3){
    allGroupData$V3<-"all"
  }
  
  result<-NULL
  for(title in unique(allGroupData$V3)){
    groupData<-allGroupData[allGroupData$V3 == title,]
    
    if(useLeastGroups){
      groupData<-groupData[which(groupData$V1 %in% samples),]
      groups<-lapply(unique(groupData$V2), function(x){
        nrow(groupData[groupData$V2==x,])
      })
      discardGroups<-NULL
      groupNames=unique(groupData$V2)
      for(i in c(1:length(groupNames))){
        sampleI<-groupData[groupData$V2==groupNames[i], "V1"]
        for(j in c(i+1:length(groupNames))){
          sampleJ<-groupData[groupData$V2==groupNames[j], "V1"]
          if(all(sampleI %in% sampleJ)){
            discardGroups<-c(discardGroups, groupNames[i])
            break
          }else if(all(sampleJ %in% sampleI)){
            discardGroups<-c(discardGroups, groupNames[j])
          }
        }
      }
      groupData<-groupData[!(groupData$V2 %in% discardGroups),]
    }
    groupData$V2<-factor(groupData$V2)
    
    res<-NULL
    gnameChanged<-FALSE
    for(sample in samples){
      stog<-groupData[groupData$V1==sample,,drop=F]
      if(nrow(stog) == 1){
        group<-stog[1,2]
      }else if(nrow(stog) > 1){
        groups<-stog$V2[order(stog$V2)]
        group<-paste(groups, collapse=":")
        gnameChanged<-TRUE
      }else{
        group<-"Unknown"
        gnameChanged<-TRUE
      }
      res<-rbind(res, data.frame(V1=sample, V2=group, V3=title))
    }
    
    if (onlySamplesInGroup) {
      #remvoe "Unknown" group
      res<-res[which(res$V2!="Unknown"),]
    }
    result<-rbind(result, res)
  }
  
  return(result)
}

gs<-getSampleInGroup(parSampleFile1, colnames(comp), onlySamplesInGroup=TRUE)
comp<-comp[,gs$V1]

rownames(gs)<-gs$V1
group<-gs[colnames(comp), "V2"]

comp_data<-comp[7:15,]
comp_data<-t(comp_data)/apply(comp_data,2,sum)

##different from any groups. 
data<-data.frame(group=group)
data$sample=DR_data(comp_data)
model1<-DirichReg(sample~group,data)
model2<-DirichReg(sample~1,data )
ares<-anova(model1,model2)

df <- data.frame(matrix(unlist(ares[1:6]), nrow=length(ares[1:6]), byrow=T))
rownames(df)<-names(ares)[1:6]
colnames(df)<-c("model_base", "model_group")
write.csv(df, file=paste0(outFile, ".anova.csv"))

plotdata<-data.frame(comp_data)
plotdata$Sample<-rownames(plotdata)
plotdata$Group<-group
mplotdata<-melt(plotdata,id.vars=c("Sample", "Group"), variable.name = "smallRNA", value.name="Proportion")

png(paste0(outFile, ".boxplot.png"), width=3000, height=2000, res=300)
g<-ggplot(mplotdata, aes(x=Group, y=Proportion, color=Group)) +   geom_boxplot() + facet_wrap(~smallRNA, scales="free_y") + theme_bw() + theme(strip.background = element_blank())
print(g)
dev.off()

png(paste0(outFile, ".heatmap.png"), width=2000, height=2000, res=300)
pheatmap(t(comp_data))
dev.off()
