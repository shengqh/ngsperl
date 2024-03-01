rm(list=ls()) 
outFile='P10473'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''

setwd('/nobackup/brown_lab/projects/20231214_10473_Methylation_hg38/MethylKitCorr/result')

### Parameter setting end ###

require(GenomicRanges)
require(tidyverse)
require(methylKit)
require(ggpubr)
require(Cairo)
library(dplyr)
library(limma)
library(minfi)
library(matrixStats)

params <- read.table(parSampleFile2, sep = "\t", header = F)
params <- params %>% column_to_rownames("V2") %>% t() %>% data.frame()
project <- params$task_name
assembly = params$org
mincov = as.numeric(params$mincov)
if("control_group" %in% names(params)){
  control_group=params$control_group
  if(control_group == ""){
    control_group=NA
  }
}else{
  control_group=NA
}

#prepare to find all specific format files
cpg.all <-read.table(parSampleFile1, sep="\t", header=F)
cpg.infile = cpg.all$V1
file.id = cpg.all$V2

meta <- read.table(parSampleFile3, sep = "\t", header=F)
colnames(meta) <- c("sample", "group")
meta<-meta[!duplicated(meta$sample),]
rownames(meta)<-meta$sample
meta<-meta[file.id,]
var = "group"

groups=unique(meta$group)
if(!is.na(control_group)){
  stopifnot(control_group %in% groups)
  groups=c(control_group, groups[groups != control_group])
}
gindex=c(1:length(groups)) - 1
names(gindex)=groups
treatment<-unlist(gindex[meta$group])

rds_file=paste0(project, ".filtered.cpg.meth.rds")
if(!file.exists(rds_file)){
  cpg.obj <- methRead(location = as.list(cpg.infile),
                      sample.id = as.list(file.id),
                      assembly = assembly,
                      treatment = treatment,
                      context = "CpG",
                      mincov = mincov)

  filtered.obj <- filterByCoverage(cpg.obj,
                                  lo.count=mincov,
                                  lo.perc=NULL,
                                  hi.count=NULL,
                                  hi.perc=99.99)
  rm(cpg.obj)

  #Merging samples
  filtered.cpg.meth <- methylKit::unite(filtered.obj, destrand=FALSE)
  rm(filtered.obj)

  #names(filtered.cpg.meth@treatment)
  saveRDS(filtered.cpg.meth, rds_file)
}else{
 filtered.cpg.meth<-readRDS(rds_file)
}

cat("PCASamples\n")
png(paste0(project, ".CpG.pca.png"), width=8, height=8, bg="white", res=300, units="in")
PCASamples(filtered.cpg.meth)
dev.off()

get_beta_value<-function(filtered.cpg.meth){
  mat = getData(filtered.cpg.meth)  
  cpg_bvalue_df = mat[, filtered.cpg.meth@numCs.index]/ (mat[,filtered.cpg.meth@numCs.index] + mat[,filtered.cpg.meth@numTs.index] )
  colnames(cpg_bvalue_df)=filtered.cpg.meth@sample.ids

  mat <- mat %>% mutate(id = paste(chr, start, end, sep = "_"))
  rownames(cpg_bvalue_df) <- mat$id

  return(cpg_bvalue_df)
}

cat("get_beta_value\n")
cpg_bvalue_df<-get_beta_value(filtered.cpg.meth)
#saveRDS(cpg_bvalue_df, paste0(project, ".cpg_bvalue_df.rds"))

hist_pdf=paste0(project, ".CpG.bvalue_hist.pdf")
if(!file.exists(hist_pdf)){
  cat("draw beta value histgram figure\n")
  pdf(hist_pdf, width=4, height=4, bg="white", onefile=TRUE)
  samples=colnames(cpg_bvalue_df)
  sample=samples[1]
  for(sample in samples){
    svalues=cpg_bvalue_df[,sample,drop=FALSE]
    g<-ggplot(svalues, aes(x = !!sym(sample))) +
      geom_histogram(bins = 100) +
      theme_bw()
    print(g)
  }
  dev.off()
}

# #https://f1000research.com/articles/6-2055/v2
# get_mds_value<-function(filtered.cpg.meth){
#   mat = getData(filtered.cpg.meth)  
#   cpg_mds_df = log2(mat[, filtered.cpg.meth@numCs.index] + 2) - log2(mat[,filtered.cpg.meth@numTs.index] + 2)
#   colnames(cpg_mds_df)=filtered.cpg.meth@sample.ids

#   id = paste0(mat$chr, "_", mat$start, "_", mat$end)
#   rownames(cpg_mds_df) <- mat$id

#   return(cpg_mds_df)
# }

# cat("get_mds_value\n")
# cpg_mds_df<-get_mds_value(filtered.cpg.meth)
# stopifnot(colnames(cpg_mds_df) == meta$sample)

# o <- order(rowVars(as.matrix(cpg_mds_df)), decreasing = TRUE)[seq_len(10000)]
# o_cpg_mds_df=cpg_mds_df[o,]

# cat("plotMDS\n")
# png(paste0(project, "_methyl_CpG_MDS_by_plotMDS.png"), width=8, height=8, bg="white", res=300, units="in")
# plotMDS(o_cpg_mds_df, top=10000, col=filtered.cpg.meth@treatment+1)
# dev.off()

#https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/mdsPlot
stopifnot(colnames(cpg_bvalue_df) == meta$sample)
png(paste0(project, ".CpG.bvalue_top10000.mdsPlot.png"), width=5, height=5, bg="white", res=300, units="in")
mdsPlot(as.matrix(cpg_bvalue_df), numPositions=10000, sampGroups=meta$group)
dev.off()

#get top 10000 most variable positions
cat("get top 10000 most variable positions\n")
bvalues_vars = matrixStats::rowVars(as.matrix(cpg_bvalue_df))
o <- order(bvalues_vars, decreasing = TRUE)[seq_len(10000)]
top10000_cpg_bvalue_df=cpg_bvalue_df[o,]

cat("cor of top10000_cpg_bvalue_df\n")
top10000_cpg_bvalue_cor <- cor(top10000_cpg_bvalue_df, method = "pearson")
saveRDS(top10000_cpg_bvalue_cor, paste0(project, ".CpG.bvalue_top10000.corr.rds"))

cat("cmdscale\n")
top10000_cpg_bvalue_corr_mds <- (1 - top10000_cpg_bvalue_cor) %>%
  cmdscale() %>%
  data.frame()
colnames(top10000_cpg_bvalue_corr_mds) <- c("Dim.1", "Dim.2")

stopifnot(rownames(top10000_cpg_bvalue_corr_mds) == colnames(cpg_bvalue_df))

top10000_cpg_bvalue_corr_mds<-top10000_cpg_bvalue_corr_mds[meta$sample,]
stopifnot(rownames(top10000_cpg_bvalue_corr_mds) == meta$sample)

top10000_cpg_bvalue_corr_mds[,var] <- meta[,var]
saveRDS(top10000_cpg_bvalue_corr_mds, paste0(project, ".CpG.bvalue_top10000.corr.MDS.rds"))

cat("plot bvalue MDS\n")
# Plot MDS
if(nrow(top10000_cpg_bvalue_corr_mds) > 20){
  cpg_label=NULL
}else{
  cpg_label=rownames(top10000_cpg_bvalue_corr_mds)
}
dms_plot <- ggscatter(top10000_cpg_bvalue_corr_mds, 
                      x = "Dim.1", 
                      y = "Dim.2",
                      label = cpg_label,
                      xlab = "MDS 1",
                      ylab = "MDS 2",
                      color = var,
                      palette = "jco",
                      font.label = 2,
                      size = 1,
                      ellipse = F,
                      ellipse.type = "norm",
                      repel = TRUE) + theme_bw() + theme(aspect.ratio=1)

ggsave(paste0(project, ".CpG.bvalue_top10000.corr.MDS.png"), dms_plot, width=5, height=4, units="in", dpi=300, bg="white")
#ggsave(paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.top10000.pdf"), dms_plot, width=4, height=3, units="in", dpi=300, bg="white")

