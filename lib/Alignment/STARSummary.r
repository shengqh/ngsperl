options(bitmapType='cairo')

library(reshape2)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputFile = args[2]

#inputFile = "/scratch/cqs/shengq1/brown/20161229_rnaseq_3685_human/star/pbs/human_3685_sample.list"
#outputFile = "/scratch/cqs/shengq1/brown/20161229_rnaseq_3685_human/star/result/human_3685.STARSummary.tsv"

cat("inputFile=", inputFile, "\n")
cat("outputFile=", outputFile, "\n")

filelist = read.table(inputFile, sep="\t", header=F, stringsAsFactors = F)

subdata = read.table(filelist$V2[1], sep="|", header=F, strip.white=T, fill=T)

final=NULL
i=1
for(i in c(1:nrow(filelist))){
  filename = filelist$V1[i]
  filelocation =filelist$V2[i]
  subdata = read.table(filelocation, sep="|", header=F, strip.white=T, fill=T, stringsAsFactors = F, row.names=1)
  colnames(subdata)= (filename)
  if(is.null(final)){
    final = subdata
  }else{
    final = cbind(final, subdata)
  }
}

write.csv(file=sub("^([^.]*).*", ".details.csv", outputFile), final)

reads=final[c("Number of input reads", "Uniquely mapped reads number", "Number of reads mapped to multiple loci", "Number of reads mapped to too many loci"),]
treads=data.frame(t(data.matrix(reads)))
write.csv(file=outputFile, treads)

colnames(treads)=c("Total", "Unique", "Multiple1", "Multiple2")
treads$Multiple=treads$Multiple1+treads$Multiple2
treads$Unmapped=treads$Total-treads$Unique-treads$Multiple
treads=treads[,c(2,5,6)]
treads$Sample=rownames(treads)

meltreads=melt(treads, id="Sample", variable.name="Read", value.name="Count")

width=max(2000, 50 * nrow(treads))
png(file=paste0(outputFile, ".png"), height=1500, width=width, res=300)
g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
  geom_bar(stat="identity", width=0.5) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
        axis.text.y = element_text(size=11, face="bold"))
print(g)
dev.off()
  
