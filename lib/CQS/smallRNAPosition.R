args = commandArgs(trailingOnly=TRUE)
positionFile = args[1]
file=tools::file_path_sans_ext(positionFile)
print(file)

fp=read.table(positionFile, sep="\t", header=T)

cols=c(rgb(139, 35, 35, 0,maxColorValue=255), 
       rgb(139, 35, 35, 255*0.2,maxColorValue=255), 
       rgb(139, 35, 35, 255*0.4,maxColorValue=255), 
       rgb(139, 35, 35,255*0.6,maxColorValue=255), 
       rgb(139, 35, 35,255*0.8,maxColorValue=255), 
       rgb(139, 35, 35,255,maxColorValue=255))

features=as.character(unique(fp$Feature))
if(length(features) > 100){
  features = features[1:100]
  fp=fp[fp$Feature %in% features,]
}
counts=lapply(features, function(x){
  idx=which(fp$Feature==x)
  unlist(fp$Count[idx])[1]
})
counts=round(unlist(counts))
strands=lapply(features, function(x){
  idx=which(fp$Feature==x)
  unlist(fp$Strand[idx])[1]
})
strands=unlist(strands)
fm=data.frame(feature=features, count=counts, strand=strands)
tnames=apply(as.matrix(fm),1,function(x){
  paste0(x[1], "(", as.numeric(x[2]),"):",x[3])
})
index = nrow(fm) + 1 - as.numeric(rownames(fm))
names(index) = features
fpf=as.character(fp$Feature)

fp$Index = index[fpf]
fp$Color=cols[round(fp$Percentage * 5)+1]

height=max(length(features) * 60, 2000)
png(paste0(positionFile, ".png"), width=8000, height=height, res=300)
plot(fp$Position, fp$Index, pch=19, cex=fp$Percentage, xlim=c(-30, 120), col=fp$Color, xlab="Position", ylab="Name and Count", yaxt="n", bty="n", xaxt="n", main=file)
axis(side=1, at=c(-10:120), las=2)
lines(x=c(0,0), y=c(0, nrow(fm)),col="blue")
text(rep(-28,nrow(fm)), index, tnames, adj=c(0,0.5))
dev.off()
