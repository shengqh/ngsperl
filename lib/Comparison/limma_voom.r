rm(list=ls()) 
outFile='LD_8511_human_cell'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/scratch/jbrown_lab/shengq2/projects/20220725_rnaseq_8511_hg38/genetable/result/LD_8511_human_cell.proteincoding.count'
parFile2=''
parFile3=''


setwd('/scratch/jbrown_lab/shengq2/projects/20220725_rnaseq_8511_hg38/limma/result')

source("countTableVisFunctions.R")

### Parameter setting end ###

options(bitmapType='cairo')

library("edgeR")
library(ggplot2)
library(patchwork)

if (grepl(".csv$",parFile1)) {
  data<-read.csv(parFile1,header=T,row.names=1,as.is=T,check.names=FALSE)
} else {
  data<-read.delim(parFile1,header=T,row.names=1,as.is=T,check.names=FALSE, sep="\t")
}

data<-data[,colnames(data) != "Feature_length"]
colClass<-sapply(data, class)
countNotNumIndex<-which(colClass!="numeric" & colClass!="integer")
if (length(countNotNumIndex)==0) {
  index<-1;
  indecies<-c()
} else {
  index<-max(countNotNumIndex)+1
  indecies<-c(1:(index-1))
}

countData<-data[,c(index:ncol(data))]
countData[is.na(countData)] <- 0
countData<-round(countData)

countAnnotations<-data[indecies,]

d0 <- DGEList(countData)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

snames <- colnames(countData) 

cultivar <- factor(gsub("_.*", "", snames))
time <- factor(gsub("_[^_]+$", "", gsub("^[^_]*_", "", snames)), levels=c("NoTreat", "TNF4hr", "TNF48hr"))
cid <- factor(gsub(".*_", "", snames))

group <- interaction(cultivar, time)
png(paste0(outFile, ".mds.png"), width=2000, height=2000, res=300)
plotMDS(d, col = as.numeric(group))
dev.off()

mm <- model.matrix(~cultivar*time)
colnames(mm)

y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
fit_coef<-coef(fit)

voomed<-as.matrix(y)

plotit <- function(pdata, gene) {
  g1<-ggplot(pdata, aes(x=cultivar, y=count)) + 
    geom_violin() + geom_boxplot(width=0.2) + 
    geom_jitter(size=1.5, width=0.1) +
    facet_wrap(~ time) + 
    stat_summary(aes(x=cultivar, y=count, group=time), fun=mean, geom="line", colour="red", size=0.8) + 
    xlab("Condition") + ylab("Voom transformed count") + ggtitle(gene) + theme_bw3()
  g2<-ggplot(pdata, aes(x=time, y=count)) + 
    geom_violin() + geom_boxplot(width=0.2) + 
    geom_jitter(size=1.5, width=0.1) +
    facet_wrap(~ cultivar) + 
    stat_summary(aes(x=time, y=count, group=cultivar), fun=mean, geom="line", colour="red", size=0.8) + 
    xlab("Time") + ylab("Voom transformed count") + ggtitle(gene) + theme_bw3()
  g<-g1+g2+plot_layout(ncol=1)
  return(g)
}

findex<-5
for(findex in c(5,6)){
  fname=gsub(":",".",colnames(fit_coef)[findex])

  tmp <- contrasts.fit(fit, coef = findex)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)

  resultAllOut<-cbind(data[rownames(top.table),], top.table)
  write.csv(resultAllOut, paste0(fname, ".limma_voom.csv"))

  sig_table<-resultAllOut[resultAllOut$adj.P.Val < 0.05 & abs(resultAllOut$logFC) > log2(1.5), ]
  write.csv(sig_table, paste0(fname, ".limma_voom.sig.csv"))

  pdf(paste0(fname, ".limma_voom.siggenes.pdf"), onefile=TRUE)
  gene_id<-rownames(sig_table)[1]
  for(gene_id in rownames(sig_table)){
    gene=sig_table[gene_id, "Feature_gene_name"]
    gdata=voomed[gene_id,]
    pdata<-data.frame("count"=gdata, "cultivar"=cultivar, "time"=time)
    g<-plotit(pdata, gene)
    print(g)
  }
  dev.off()
}
