#readFile = "CR_Y_TBX5_peaks.broadPeak.bed.reads" 
#singlePdf = 0 
#inputFile = "CR_Y_TBX5_peaks.broadPeak.bed.depth" 
#outputFile = "" 
#facet<-0
#drawLine<-1

library("reshape2")
library("ggplot2")

data<-read.table(inputFile, sep="\t", header=T, stringsAsFactors = F)

if(exists("readFile")){
  sampleReads<-read.table(readFile, sep="\t", header=T, row.names=1, as.is = T)
  totalReads<-sampleReads[, 1]
  names(totalReads)<-rownames(sampleReads)
  for(sample in rownames(sampleReads)){
    data[,sample] = data[,sample] * 1000000 / totalReads[sample]
  }
}

if(exists("cnvrFile")){
  cnvr<-read.table(cnvrFile, sep="\t", header=T, stringsAsFactors = F, row.names=4)
  refs<-rownames(sampleReads)[! rownames(sampleReads) %in% colnames(cnvr)]
  
  colors<-c("green", "darkblue", "lightblue", "black", colorRampPalette(c("yellow", "red"))(11))
  names(colors)<-c("REF","CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN16", "CN32", "CN64")
  
  #no_sig<-c("CN1","CN2","CN3","REF")
  no_sig<-c("CN2","REF")
}

files<-unique(data$File)

if(singlePdf){
  pdf(outputFile, onefile = T)
}

x<-files[2]
for(x in files){
  cat(x, "\n")
  
  if(exists("cnvrFile")){
    tmpcnv<-cnvr[x, c(4:ncol(cnvr))]
    tmpcnv[,refs] <- "REF"
    tmpcnv<-t(tmpcnv)
    
    curcnv<-as.character(tmpcnv[,1])
    names(curcnv) <- row.names(tmpcnv)
    
    if(length(curcnv[! (curcnv %in% no_sig)]) == 0){
      next
    }
  }
  
  curdata<-data[data$File==x,]
  
  title<-paste0(x, " (", curdata$Chr[1], ":", min(curdata$Position),"-",max(curdata$Position),")")
  
  mdata<-melt(curdata, id=c("Chr", "Position", "File"))
  colnames(mdata)<-c("Chr", "Position", "File", "Sample", "Depth")
  
  if(facet){
    height=max(2000, 400+300 * length(unique(mdata$Sample)))
    if(exists("cnvrFile")){
      mdata$Color<-as.character(curcnv[as.character(mdata$Sample)])
      g<-ggplot(mdata, aes(x=Position, y=Depth))
      if(drawLine){
        g <- g + geom_line(aes(color = Color), size=0.8)
      }else{
        g <- g + geom_point(aes(color = Color), size=0.8)
      }
      g<-g + scale_colour_manual(name="CNV", values = colors)
    }else{
      g<-ggplot(mdata, aes(x=Position, y=Depth))
      if(drawLine){
        g<-g+geom_line(aes(color = Sample), size=0.8, show.legend = F)
      }else{
        g<-g+geom_point(aes(color = Sample), size=0.8, show.legend = F)
      }
    }
    g <- g + xlab(unique(data$chr)) + 
      ylab("Reads per million total reads") +
      ggtitle(x) +
      facet_wrap( ~ Sample, ncol=1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  else{
    height=2000
    if(exists("cnvrFile")){
      mdata$Color<-as.character(curcnv[as.character(mdata$Sample)])
      g<-ggplot(mdata, aes(x=Position, y=Depth, group=Sample))
      if(drawLine){
        g <- g + geom_line(aes(color = Color), size=0.8)
      }else{
        g <- g + geom_point(aes(color = Color), size=0.8)
      }
      g<-g+scale_colour_manual(name="CNV", values = colors)
    }else{
      g<-ggplot(mdata, aes(x=Position, y=Depth, group=Sample))
      if(drawLine){
        g <- g + geom_line(aes(color = Sample), size=0.8, show.legend = T)
      }else{
        g <- g + geom_point(aes(color = Sample), size=0.8, show.legend = T)
      }
    }
    
    g <- g + xlab(data$chr[1]) + 
      ylab("Reads per million total reads") +
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  if(singlePdf){
    print(g)
  }else{
    png(paste0(x, ".png"), width=2000, height=height, res=300)
    print(g)
    dev.off()
  }
}

if(singlePdf){
  dev.off()
}
