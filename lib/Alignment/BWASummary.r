rm(list=ls()) 
outFile='human_exomeseq'
parSampleFile1=''
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile6='fileList6.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20210616_human_exomeseq/bwa_g4_refine_summary/result')

### Parameter setting end ###

source("AlignmentUtils.r")
options(bitmapType='cairo')

library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)

#parSampleFile1 = "/scratch/cqs/shengq2/dnaseq/20201014_liuqi_gene_panel_dbsnp150_bigbed/bwa_summary/result/Adenoma__fileList1.list"
#outFile = "/scratch/cqs/shengq2/dnaseq/20201014_liuqi_gene_panel_dbsnp150_bigbed/bwa_summary/result/Adenoma"

#source("AlignmentUtils.r")

if(!exists("rg_name_regex")){
  rg_name_regex=NA
}

if (parSampleFile2 != ""){
  draw_chromosome_count(listFile=parSampleFile2, 
                        outFilePrefix=outFile,
                        rg_name_regex=rg_name_regex, 
                        remove_chrM_genes=FALSE,
                        draw_read_png=parSampleFile1 == "")
}

if (parSampleFile1 != ""){
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

  if (!is.na(rg_name_regex)) {
    sample_names<-str_match(colnames(final), rg_name_regex)[,2]
    unique_sample_names<-unique(sample_names)

    for (sample_name in unique_sample_names){
      columns = which(sample_names == sample_name)
      col_eqn <- paste(colnames(final)[columns], collapse = " + ")
      newdata = final %>% mutate(combined = eval(parse(text = col_eqn)))
      final[, sample_name] = newdata$combined
    }
    final<-final[,unique_sample_names,drop=F]
  }
  write.csv(file=paste0(outFile, ".BWASummary.details.csv"), final)

  reads=final[c("UnmappedFragments", "UniquelyMappedFragments", "MultipleMappedFragments"),]
  treads=data.frame(t(data.matrix(reads)))
  write.csv(file=paste0(outFile, ".BWASummary.csv"), treads)

  colnames(treads)=c("Unmapped", "Unique", "Multiple")
  treads$Sample=rownames(treads)

  meltreads=melt(treads, id="Sample", variable.name="Read", value.name="Count")
  meltreads$Read<-factor(meltreads$Read, levels=c("Multiple", "Unmapped", "Unique"))

  width=max(2000, 50 * nrow(treads))

  draw_reads<-function(meltreads, filename, width){
    png(file=filename, height=1500, width=width, res=300)
    g=ggplot(meltreads, aes(x=Sample, y=Count, fill=Read)) + 
      geom_bar(stat="identity", width=0.5) +
      theme(axis.text.x = element_text(angle=90, vjust=0.5, size=11, hjust=0, face="bold"),
            axis.text.y = element_text(size=11, face="bold")) + theme_bw() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("") + ylab("Reads")
    print(g)
    dev.off()
  }
  draw_reads(meltreads, paste0(outFile, ".reads.png"), width)

  treads$Total<-apply(treads[,c(1:3)], 1, sum)
  treads$UnmappedPerc<-treads$Unmapped / treads$Total
  treads<-treads[order(treads$UnmappedPerc, decreasing=T),]

  meltreads$Sample<-factor(meltreads$Sample, levels=treads$Sample)

  draw_reads(meltreads, paste0(outFile, ".reads.sorted.png"), width)
}
