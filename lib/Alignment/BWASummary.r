options(bitmapType='cairo')

library(reshape2)
library(ggplot2)

#parSampleFile1 = "/scratch/cqs/shengq2/dnaseq/20201014_liuqi_gene_panel_dbsnp150_bigbed/bwa_summary/result/Adenoma__fileList1.list"
#outFile = "/scratch/cqs/shengq2/dnaseq/20201014_liuqi_gene_panel_dbsnp150_bigbed/bwa_summary/result/Adenoma"

filelist = read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)

final=NULL
i=1
for(i in c(1:nrow(filelist))){
  filename = filelist$V2[i]
  filelocation =filelist$V1[i]
  subdata = read.table(filelocation, sep="\t", header=T, strip.white=T, fill=T, stringsAsFactors = F, row.names=1)
  colnames(subdata)= (filename)
  if(is.null(final)){
    final = subdata
  }else{
    final = cbind(final, subdata)
  }
}

write.csv(file=paste0(outFile, ".BWASummary.details.csv"), final)

reads=final[c("UnmappedFragments", "UniquelyMappedFragments", "MultipleMappedFragments"),]
treads=data.frame(t(data.matrix(reads)))
write.csv(file=paste0(outFile, ".BWASummary.csv"), treads)

colnames(treads)=c("Unmapped", "Unique", "Multiple")
treads$Sample=rownames(treads)

meltreads=melt(treads, id="Sample", variable.name="Read", value.name="Count")

width=max(2000, 50 * nrow(treads))
png(file=paste0(outFile, ".BWASummary.png"), height=1500, width=width, res=300)
g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
  geom_bar(stat="identity", width=0.5) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
        axis.text.y = element_text(size=11, face="bold"))
print(g)
dev.off()

treads$Total<-apply(treads[,c(1:3)], 1, sum)
treads$UnmappedPerc<-treads$Unmapped / treads$Total
treads<-treads[order(treads$UnmappedPerc, decreasing=T),]

meltreads$Sample<-factor(meltreads$Sample, levels=treads$Sample)

png(file=paste0(outFile, ".BWASummary.sorted.png"), height=1500, width=width, res=300)
g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
  geom_bar(stat="identity", width=0.5) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
        axis.text.y = element_text(size=11, face="bold"))
print(g)

  
