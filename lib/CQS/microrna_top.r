#predefine_start

outputdir<-"Z:/Shared/Labs/Vickers Lab/ShilinZhao/20151017_3018-KCV-45-46/bowtie1_genome_1mm_NTA_smallRNA_table/result"
inputfile<-"smallRNA_1mm_3018-KCV-45-46.miRNA.count"
outputfile<-"smallRNA_1mm_3018-KCV-45-46.miRNA.count.png"
replaceSampleName<-0

#predefine_end

options(bitmapType='cairo')

setwd(outputdir)

library(ggplot2)
library(reshape2)

data<-read.table(inputfile, header=T, sep="\t", stringsAsFactor=F, row.names=1)
isDataNumeric = unlist(lapply(data[1,], function(x){is.numeric(x)}))
if (any(isDataNumeric)) {
  index = 1
  while(!all(isDataNumeric[index:ncol(data)])){
    index = index + 1
  }
} else {
  cat("Error: No numeric data found for ", inputfile, "\n")
  quit(save="no")
}

if(index > 1){
  indecies<-c(1:(index-1))
}else{
  indecies<-c()
}
countData<-data[,c(index:ncol(data)),drop=F]
rownames(countData)<-gsub(";.*", "", rownames(countData))
if(replaceSampleName){
  colnames(countData)<-paste0("S", c(1:ncol(countData)))
}

totalCount<-apply(countData,2,sum)

percData<-data.frame(prop.table(as.matrix(countData), margin=2)*100)

top3indecies<-unique(melt(data.frame(apply(percData, 2, function(x){
  unlist(order(x, decreasing=TRUE)[1:3])
})))$value)

top3indecies<-sort.int(top3indecies)
perctable<-percData[top3indecies,,drop=F]
perctable$Feature<-factor(rownames(perctable), levels=rownames(perctable))
melttable<-melt(perctable, id="Feature")
colnames(melttable)<-c("Feature", "Sample", "Percentage")

png(filename = outputfile, height=1000, width=max(2000, ncol(countData) * 100), res=300)
g<-ggplot(melttable, aes(x=Sample, y=Percentage, fill=Feature))+
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=1, face="bold"),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))
print(g)
dev.off()
