require("ggplot2")
require("data.table")
require("stringr")

draw_chromosome_count<-function(listFile, outFilePrefix, rg_name_regex=NA, remove_chrM_genes=FALSE) {
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
  all_chroms = c()
  i=1
  for(i in c(1:nrow(filelist))){
    filename = filelist$V2[i]
    filelocation =filelist$V1[i]

    if(!file.exists(filelocation)){
      missing = c(missing, filelocation)
      next
    }

    cat("reading", filelocation, "\n")
    subdata = fread(filelocation, data.table=FALSE)
    if(colnames(subdata)[1] == "Geneid") { 
      #featureCounts count file, we need to aggregate it to chromosome level
      subdata = subdata[,c(2,7)]
      subdata$Chr = gsub(";.+","",subdata$Chr)
      colnames(subdata)[2]="Reads"
      subdata = aggregate(subdata$Reads, by=list(subdata$Chr), FUN=sum)
      colnames(subdata)=c("Chrom", "Reads")
    }else{
      subdata = read.table(filelocation, sep="\t", header=F, stringsAsFactors = F)
      subdata=subdata[,c(1,3)]
      colnames(subdata)<-c("Chrom", "Reads")
    }
    subdata=subdata[subdata$Chrom != '*',]
    subdata$Sample=filename

    chrdata=subdata[str_length(subdata$Chrom) < 6,]
    mediancount<-median(chrdata$Reads)
    sdata<-subdata[str_length(subdata$Chrom) < 6 | subdata$Reads > mediancount,]
    all_chroms<-c(all_chroms, sdata$Chrom)
    final=rbind(final, subdata)
  }
  all_chroms<-sort(unique(all_chroms))
  if(remove_chrM_genes){
    all_chroms=all_chroms[all_chroms != "chrM"]
  }

  final=final[final$Chrom %in% all_chroms,]

  if(length(missing) > 0){
    writeLines(missing, paste0(outFilePrefix, ".chromosome.missing"))
    chromosomeFilePrefix = paste0(outFilePrefix, ".chromosome.valid")
  }else{
    chromosomeFilePrefix = paste0(outFilePrefix, ".chromosome")
  }

  is_bam_scattered<-!is.na(rg_name_regex)

  if(is_bam_scattered){
    write.csv(file=paste0(chromosomeFilePrefix, ".scattered.csv"), final, row.names=F)
  }else{
    write.csv(file=paste0(chromosomeFilePrefix, ".csv"), final, row.names=F)
  }

  chroms=paste0("chr", c(1:22,'X','Y','M', 'MT'))
  if(!any(all_chroms %in% chroms)){
    chroms=paste0("", c(1:22,'X','Y','M', 'MT'))
  }
  chroms=chroms[chroms %in% all_chroms]
  chroms<-c(chroms, all_chroms[!(all_chroms %in% chroms)])

  fc<-reshape2::dcast(final, "Sample ~ Chrom", value.var="Reads", fill=0)
  final<-reshape2::melt(fc, id="Sample")
  colnames(final)<-c("Sample", "Chrom", "Reads")
  final$Chrom=factor(final$Chrom, levels=chroms)

  final$NoRead=final$Reads==0

  noreads<-final[final$NoRead,]
  noreads<-noreads[!duplicated(noreads$Sample),]
  if(is_bam_scattered){
    write.table(file=paste0(chromosomeFilePrefix, ".scattered.noread.txt"), noreads, sep="\t", quote=F, row.names=F)
  }else{
    write.table(file=paste0(chromosomeFilePrefix, ".noread.txt"), noreads, sep="\t", quote=F, row.names=F)
  }

  colors=c("black","red")
  names(colors)=c("FALSE","TRUE")

  draw_figure<-function(final, filename){
    if(any(final$NoRead)){
      g<-ggplot(final, aes(x=Chrom, y=Sample)) + 
        geom_point(aes(size=Reads, color=NoRead)) + theme_classic() + 
        scale_color_manual(values=colors) +
        theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0),
              axis.title = element_blank())
    }else{
      g<-ggplot(final, aes(x=Chrom, y=Sample)) + 
        geom_point(aes(size=Reads)) + theme_classic() + 
        theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0),
              axis.title = element_blank())
    }

    fontsize_inch = GeomLabel$default_aes$size * 0.0393701
    height=max(4, fontsize_inch * length(unique(final$Sample)) + 1)
    ggsave(filename, g, width=10, height=height, units="in", dpi=300, bg="white", limitsize = FALSE)
  }

  if(!is.na(rg_name_regex)){
    final2<-final[,c("Sample", "Chrom", "Reads")]
    final2$Sample<-str_match(final2$Sample, rg_name_regex)[,2]
    final3 <- aggregate(Reads ~ Sample + Chrom, data = final2, FUN = sum, na.rm = TRUE)
    final3$NoRead=final3$Reads==0
    draw_figure(final3, paste0(chromosomeFilePrefix, ".png"))
    write.csv(final3, paste0(chromosomeFilePrefix, ".csv"))
  }else{
    draw_figure(final, paste0(chromosomeFilePrefix, ".png"))
  }
}
