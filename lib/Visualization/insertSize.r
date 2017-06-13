options(bitmapType='cairo')
library(Rsamtools)
library(ggplot2)

filelist<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactor=F)
fileDensities<-apply(filelist, 1, function(x){
  res<-scanBam(x[1], index=x[1], param=ScanBamParam(what=c("isize")))
  aa<-res[[1]]
  bb<-as.numeric(aa[[1]])
  cc<-abs(bb)
  cc<-cc[!is.na(cc) & cc < 1000]
  dens<-density(cc, from=0, to=1000, n=250)
  dd<-data.frame(File=x[2], InsertSize=dens$x, Density=dens$y )
  return(dd)
})
all<-do.call(rbind, fileDensities)
png(file=outFile, width=2000, height=2000, res=300)
g<-ggplot(all, aes(InsertSize, Density, group=File, color=File)) + geom_line() + theme_bw() + scale_x_continuous(breaks=seq(0,1000,100))
print(g)
dev.off()
