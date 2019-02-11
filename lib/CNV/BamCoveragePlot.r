library(ggplot2)
library(reshape2)

dat<-read.table(covfile, sep="\t", header=F, stringsAsFactors = F)

sampleIndex<-4
for (sampleIndex in c(4:ncol(dat))){
  sampleData<-dat[,c(1,2,sampleIndex)]
  sampleName<-sample_names[sampleIndex-3]
  colnames(sampleData)<-c("Chromosome", "Position", "Reads")
  sampleData<-sampleData[sampleData$Chromosome != "Y",]
  sampleData$Chromosome[sampleData$Chromosome == "X"]<-"23"
  sampleData$Chromosome<-as.numeric(sampleData$Chromosome)
  
  finalData<-sampleData[FALSE,]
  
  chr<-sampleData$Chromosome[1]
  for(chr in unique(sampleData$Chromosome)){
    chrData<-sampleData[sampleData$Chromosome==chr,]
    csumz<-which(chrData$Reads>0)
    if(length(csumz) == 0){
      next
    }
    
    mincsumz<-min(csumz)
    maxcsumz<-max(csumz)

    chrData<-chrData[c(mincsumz:maxcsumz),]
    finalData<-rbind(finalData, chrData)
  }
  finalData$Log2Reads<-log2(finalData$Reads+1)
  
  pdf(file=paste0(covfile, "_", sampleName, ".log2.pdf"), width=5, height=20)
  g<-ggplot(finalData, aes(x=Position, y=Log2Reads)) + 
    geom_point(size=0.5) + 
    facet_grid(Chromosome~.) +
    ylab("log2(Reads)") +
    theme_classic()
  print(g)
  dev.off()
}
