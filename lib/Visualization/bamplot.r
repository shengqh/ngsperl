library(reshape2)
library(ggplot2)

args <- commandArgs(TRUE)

print(args)

summaryfile = args[1]
width=as.numeric(args[2])
height=as.numeric(args[3])

#summary=read.table("KCV_3018_77_78_79/KCV_3018_77_78_79_summary.txt", header=T, sep="\t")
summary=read.table(summaryfile, header=T, sep="\t", stringsAsFactor=F)
summary$title = paste0(summary$CHROM, ":", summary$START, "-", summary$END)

idx=1
for (idx in c(1:nrow(summary))){
  plotfile=summary$PLOT_TABLE[idx]
  id=summary$ID[idx]
  ptitle=summary$title[idx]
  
  plotdata<-read.table(plotfile, header=T, sep="\t")

  pdata<-plotdata[,c(2, 3, 8:ncol(plotdata))]
  mdata<-melt(pdata, id=c("GENE_ID", "NAME"), value.name="Reads")
  mdata$variable<-gsub("bin_", "", mdata$variable)
  mdata$variable<-as.numeric(mdata$variable)
  
  pdf(paste0(id, "_", gsub("[:-]", "_", ptitle), ".pdf"), width=width, height=height)
  g=ggplot(mdata, aes(x=variable, y=Reads)) + 
    geom_line() + 
    facet_wrap(~NAME, strip.position="bottom") + 
    xlab(ptitle) + 
    ylim(0, NA) +
    theme_classic() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background=element_blank())
  print(g)
  dev.off()
}

