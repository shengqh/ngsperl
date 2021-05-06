require("ggplot2")
require("data.table")
require("stringr")

draw_chromosome_count<-function(listFile, outFilePrefix) {
  filelist = read.table(listFile, sep="\t", header=F, stringsAsFactors = F)

  final=NULL
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]
    subdata = read.table(filelocation, sep="\t", header=F, stringsAsFactors = F)
    subdata=subdata[,c(1,3)]
    colnames(subdata)<-c("Chrom", "Reads")
    subdata=subdata[str_length(subdata$Chrom) < 6,]
    subdata=subdata[subdata$Chrom != '*',]
    subdata$Sample=filename
    final=rbind(final, subdata )
  }
  write.csv(file=paste0(outFilePrefix, ".chromosome.csv"), final, row.names=F)

  chroms=paste0("chr", c(1:22,'X','Y','M', 'MT'))
  if(!any(final$Chrom %in% chroms)){
    chroms=paste0("", c(1:22,'X','Y','M', 'MT'))
  }
  chroms=chroms[chroms %in% final$Chrom]
  final$Chrom=factor(final$Chrom, levels=chroms)
  final$NoRead=final$Reads==0

  colors=c("black","red")
  names(colors)=c("FALSE","TRUE")

  height=max(1000, 60 * length(unique(final$Sample)))
  png(file=paste0(outFilePrefix, ".chromosome.png"), height=height, width=3000, res=300)
  g<-ggplot(final, aes(x=Chrom, y=Sample)) + 
    geom_point(aes(size=Reads, color=NoRead)) + theme_classic() + 
    scale_color_manual(values=colors) +
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0),
          axis.title = element_blank())
  print(g)
  dev.off()
}
