require("ggplot2")
require("data.table")
require("stringr")

draw_chromosome_count<-function(listFile, outFilePrefix) {
  filelist = read.table(listFile, sep="\t", header=F, stringsAsFactors = F)

  missing = c()
  missing_file=paste0(outFilePrefix, ".chromosome.missing")
  valid_csv=paste0(outFilePrefix, ".chromosome.valid.csv")
  valid_png=paste0(outFilePrefix, ".chromosome.valid.png")
  if(file.exists(missing_file)){
    file.remove(missing_file)
  }
  if(file.exists(valid_csv)){
    file.remove(valid_csv);
  }
  if(file.exists(valid_png)){
    file.remove(valid_png);
  }

  filelist$size=file.info(filelist$V1)$size

  missing = filelist$V1[is.na(filelist$size) | filelist$size==0]
  filelist=filelist[!(filelist$V1 %in% missing),]

  final=NULL
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]

    if(!file.exists(filelocation)){
      missing = c(missing, filelocation)
      next
    }

    subdata = read.table(filelocation, sep="\t", header=F, stringsAsFactors = F)
    subdata=subdata[,c(1,3)]
    colnames(subdata)<-c("Chrom", "Reads")
    subdata=subdata[str_length(subdata$Chrom) < 6,]
    subdata=subdata[subdata$Chrom != '*',]
    subdata$Sample=filename
    final=rbind(final, subdata )
  }

  if(length(missing) > 0){
    writeLines(missing, paste0(outFilePrefix, ".chromosome.missing"))
    chromosomeFilePrefix = paste0(outFilePrefix, ".chromosome.valid")
  }else{
    chromosomeFilePrefix = paste0(outFilePrefix, ".chromosome")
  }

  write.csv(file=paste0(chromosomeFilePrefix, ".chromosome.csv"), final, row.names=F)

  chroms=paste0("chr", c(1:22,'X','Y','M', 'MT'))
  if(!any(final$Chrom %in% chroms)){
    chroms=paste0("", c(1:22,'X','Y','M', 'MT'))
  }
  chroms=chroms[chroms %in% final$Chrom]

  fc<-reshape2::dcast(final, "Sample ~ Chrom", value.var="Reads", fill=0)
  final<-melt(fc, id="Sample")
  colnames(final)<-c("Sample", "Chrom", "Reads")
  final$Chrom=factor(final$Chrom, levels=chroms)

  final$NoRead=final$Reads==0

  noreads<-final[final$NoRead,]
  noreads<-noreads[!duplicated(noreads$Sample),]
  write.table(file=paste0(chromosomeFilePrefix, ".noread.txt"), noreads, sep="\t", quote=F, row.names=F)

  colors=c("black","red")
  names(colors)=c("FALSE","TRUE")

  height=max(1000, 60 * length(unique(final$Sample)))
  png(file=paste0(chromosomeFilePrefix, ".png"), height=height, width=3000, res=300)
  g<-ggplot(final, aes(x=Chrom, y=Sample)) + 
    geom_point(aes(size=Reads, color=NoRead)) + theme_classic() + 
    scale_color_manual(values=colors) +
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0),
          axis.title = element_blank())
  print(g)
  dev.off()
}
