options(bitmapType='cairo')

library(reshape2)
library(ggplot2)

filelist = read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

finalFile = paste0(outFile, ".csv")
final=NULL
i=1
for(i in c(1:nrow(filelist))){
  filename = filelist$V2[i]
  filelocation = filelist$V1[i]
  if(file.exists(filelocation)){
    subdata = read.table(filelocation, sep="\t", header=T, stringsAsFactors = F)
    subdata$Sample = filename
    final<-rbind(final, subdata)
  }else{
    finalFile = paste0(outFile, ".tmp.csv")
  }
}

write.csv(final, file=finalFile, row.names=F)

height=length(unique(final$Sample)) * 100
#cat(height)
png(file=paste0(finalFile, ".png"), width=9000, height=height, res=300)
g<-ggplot(final, aes(x=Sample, y=Count)) + 
  geom_bar(aes(fill=Category), stat = "identity") + 
  coord_flip() +
  facet_grid(rows=~Chromosome, scales="free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(g)
dev.off()
