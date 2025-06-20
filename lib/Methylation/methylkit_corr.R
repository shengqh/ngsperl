rm(list=ls()) 
outFile='P13303'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20250605_13303_DNAMethyl_mm10/MethylKitCorr/result')

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
mds_legendPos = params$mds_legendPos

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

stopifnot(colnames(cpg_bvalue_df) == meta$sample)

#https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/mdsPlot
#Use the mdsPlot function from minfi package to draw MDS plot by Euclidean distance.

png(paste0(project, ".euclidean_distance.all.MDS.png"), width=5, height=5, bg="white", res=300, units="in")
mdsPlot(as.matrix(cpg_bvalue_df), numPositions=nrow(cpg_bvalue_df), sampGroups=meta$group, pch=16, legendPos=mds_legendPos, legendNCol=1)
dev.off()

png(paste0(project, ".euclidean_distance.top10000.MDS.png"), width=5, height=5, bg="white", res=300, units="in")
mdsPlot(as.matrix(cpg_bvalue_df), numPositions=10000, sampGroups=meta$group, pch=16, legendPos=mds_legendPos, legendNCol=1)
dev.off()

draw_corr_mds_plot <- function(cpg_bvalue_df, groups, output_file=NULL){
  cat("cor ...\n")
  cpg_bvalue_cor <- cor(cpg_bvalue_df, method = "pearson")

  cat("cmdscale ...\n")
  cpg_bvalue_corr_mds <- (1 - cpg_bvalue_cor) %>%
    cmdscale() %>%
    data.frame()
  colnames(cpg_bvalue_corr_mds) <- c("Dim.1", "Dim.2")

  cpg_bvalue_corr_mds<-cpg_bvalue_corr_mds[meta$sample,]
  cpg_bvalue_corr_mds$groups = groups

  cat("plot bvalue MDS ...\n")
  # Plot MDS
  if(nrow(cpg_bvalue_corr_mds) > 20){
    cpg_label=NULL
  }else{
    cpg_label=rownames(cpg_bvalue_corr_mds)
  }
  dms_plot <- ggscatter(cpg_bvalue_corr_mds, 
                        x = "Dim.1", 
                        y = "Dim.2",
                        label = cpg_label,
                        xlab = "MDS 1",
                        ylab = "MDS 2",
                        color = "groups",
                        palette = "jco",
                        font.label = 6,
                        size = 2,
                        ellipse = F,
                        ellipse.type = "norm",
                        repel = TRUE) + theme_bw() + theme(aspect.ratio=1)

  if(!is.null(output_file)){
    ggsave(output_file, dms_plot, width=5, height=4, units="in", dpi=300, bg="white")
  }  
  return(dms_plot)
}

cpg_bvalue_df = cpg_bvalue_df[,meta$sample]

cat("draw MDS plot for all CPGs ...\n")
saved = draw_corr_mds_plot(cpg_bvalue_df, meta$group, paste0(project, ".pearson_corr.all.MDS.png"))

cat("get top 10000 most variable positions ...\n")
bvalues_vars = matrixStats::rowVars(as.matrix(cpg_bvalue_df))
o <- order(bvalues_vars, decreasing = TRUE)[seq_len(10000)]
top10000_cpg_bvalue_df=cpg_bvalue_df[o,]

cat("draw MDS plot for top 10000 CPGs ...\n")
saved = draw_corr_mds_plot(top10000_cpg_bvalue_df, meta$group, paste0(project, ".pearson_corr.top10000.MDS.png"))
