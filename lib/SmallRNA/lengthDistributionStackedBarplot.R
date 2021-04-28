
library(ggplot2)
library(data.table)

fileList1 <- read.table("fileList1.txt", header=FALSE, sep="\t", stringsAsFactors = F, row.names=2)
RNA_class <- rownames(fileList1)[!(rownames(fileList1) %in% c("category"))]

category <- read.csv(fileList1["category",1], header=T, row.names=1, check.names=F)
totalReads<-unlist(category["TotalReads",])

final<-NULL
for(rna in RNA_class){
  print(rna)
  read.count<-fread(fileList1[rna,1],data.table=FALSE,check.names=TRUE)

  if(is.numeric(read.count[,1])){
    read.count.length=read.count[,1]
  }else{
    read.count.length=nchar(read.count[,1])
  }

  if("TopFeature" %in% colnames(read.count)){
    read.count<-read.count[,c(-1,-2,-3)]
  }else{
    read.count<-read.count[,-1]
  }

  read.count.agg<-aggregate(read.count, by=list(read.count.length),sum)
  row.names(read.count.agg)<-read.count.agg[,1]
  read.count.agg<-read.count.agg[,colnames(category)]
  read.count.agg.rpm<-t(t(read.count.agg) / totalReads * 1e+6)
  
  mdata<-melt(read.count.agg.rpm)
  colnames(mdata)<-c("Length", "Sample", "RPM")
  mdata$Category<-rna
  if(is.null(final)){
    final<-mdata
  }else{
    final<-rbind(final, mdata)
  }
}

fastq_length<-final[final$Category=="fastq_len",]
final<-final[final$Category!="fastq_len",]

RNA_class<-RNA_class[RNA_class != "fastq_len"]
final$Category<-factor(final$Category, levels=(RNA_class))

allcolors<-c("blue","green","red", "brown","purple", "yellow", "black")[1:length(unique(final$Category))]

textTitle<-element_text(face= "bold", color = "black", size=22, hjust=0.5)
text20Bold<-element_text(face= "bold", color = "black", size=20)
text20<-element_text(color = "black", size=20)

pdf(file=paste0(outFile,".pdf"), onefile=T, width=8, height=7)

for(sample in unique(fastq_length$Sample)){
  fq<-fastq_length[fastq_length$Sample==sample,]
  fs<-final[final$Sample==sample,]
  g<-ggplot()+ 
    geom_area(data=fq, aes(x=Length,y=RPM),fill="gray85",color="black")+
    geom_bar(data=fs, aes(x=Length,y=RPM, fill=Category),
             stat="identity", width=0.8, color=NA,size=0.73)+
    scale_fill_manual(values=allcolors)+
    theme_classic() + 
    labs(x= "Read Length", 
         y="Total Reads Per Million",
         title=sample)+
    theme(plot.title = textTitle,
          axis.title = text20Bold,
          axis.text = text20,
          axis.line = element_line(colour = "gray75", size =0.73, linetype = "solid"),
          axis.ticks = element_line(size=0.73),axis.ticks.length=unit(0.3,"cm"),
          legend.text = text20Bold,
          legend.title=element_blank())
  print(g)  
}

dev.off()

groupFileList<-parSampleFile2

cellWidth=1500
scales="free"

sampleCount<-length(unique(fastq_length$Sample))
facetColCount = ceiling(sqrt(sampleCount))
facetRowCount = ceiling( sampleCount * 1.0 / facetColCount)
width=max(2000, cellWidth * facetColCount + 300)
height=max(2000, cellWidth * facetRowCount)

png(file=paste0(outFile, ".png"), width=width, height=height, res=300)

g<-ggplot()+ 
  geom_area(data=fastq_length, aes(x=Length,y=RPM),fill="gray85",color="black")+
  geom_bar(data=final, aes(x=Length,y=RPM, fill=Category),
           stat="identity", width=0.8, color=NA,size=0.73)+
  scale_fill_manual(values=allcolors)+
  theme_classic() + 
  labs(x = "Read Length", 
       y = "Total Reads Per Million")+
  theme(plot.title = textTitle,
        axis.title = text20Bold,
        axis.text = text20,
        axis.line = element_line(colour = "gray75", size=0.73, linetype = "solid"),
        axis.ticks = element_line(size=0.73),axis.ticks.length=unit(0.3,"cm"),
        strip.background = element_blank(),
        strip.text = text20Bold,
        legend.text = text20Bold,
        legend.title = element_blank()) +
        facet_wrap(~Sample, ncol=facetColCount, scales=scales)
print(g)

dev.off()
