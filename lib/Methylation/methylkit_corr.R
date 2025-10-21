rm(list=ls()) 
outFile='P13927'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20251006_13927_DNAMethyl_mm10/MethylKitCorr/result')

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
library(ggrepel)

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

if("pipeline" %in% names(params)){
  pipeline=params$pipeline
  if(pipeline == ""){
    pipeline="amp"
  }
}else{
  pipeline="amp"
}
cat("Using data from pipeline:", pipeline, "\n")

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
                      mincov = mincov,
                      pipeline = pipeline)

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
stopifnot(colnames(cpg_bvalue_df) == meta$sample)
#saveRDS(cpg_bvalue_df, paste0(project, ".cpg_bvalue_df.rds"))

# remove rows containing NA values, they might be introduced at unite step
cpg_bvalue_df = cpg_bvalue_df[ rowSums(is.na(cpg_bvalue_df))==0, ] 

cpg_bvalue_df_m = reshape2::melt(cpg_bvalue_df, value.name = "bvalue")
hist_ncol=ceiling(sqrt(ncol(cpg_bvalue_df)))
hist_nrow=ceiling(ncol(cpg_bvalue_df)/hist_ncol)

g<-ggplot(cpg_bvalue_df_m, aes(x = bvalue)) +
  geom_histogram(bins = 100) +
  facet_wrap(~variable, ncol=hist_ncol) +
  theme_bw() +
  xlab("Beta value") +
  ylab("CpG sites") +
  theme(strip.background=element_rect(fill="white"))

hist_png=paste0(project, ".CpG.bvalue_hist.png")
ggsave(hist_png, g, width=hist_ncol*2, height=hist_nrow*1.5, units="in", dpi=300, bg="white")

sds=rowSds(as.matrix(cpg_bvalue_df))

cat("get top 10000 most variable positions ...\n")
o <- order(sds, decreasing = TRUE)[seq_len(10000)]
top10000_cpg_bvalue_df=cpg_bvalue_df[o,]

cat("keep CpG with non-zero standard deviation ...\n")
cpg_bvalue_df=cpg_bvalue_df[sds>0,]

if(0){
  cur_bvalue_df=cpg_bvalue_df
  groups=meta$group
  output_file=paste0(project, ".CpG.all.PCA.png")
  title="CpG methylation PCA Analysis"
}

draw_pca_plot <- function(cur_bvalue_df, groups=NULL, output_file=NULL, title="CpG methylation PCA Analysis"){
  cat("prcomp ...\n")

  cpg_pr = prcomp(t(cur_bvalue_df),scale=TRUE,center=TRUE)
  supca = summary(cpg_pr)$importance

  fit = as.data.frame(scale(cpg_pr$x))
  pcalabs=paste0(colnames(fit), "(", round(supca[2,] * 100), "%)")

  fit = fit |> tibble::rownames_to_column("Sample") 

  hasGroup = !is.null(groups)
  if(hasGroup){
    fit$groups = groups
  }else{
    fit$groups = "All"
  }

  if(nrow(fit) > 20){
    cpg_label=NULL
  }else{
    cpg_label=fit$Sample
  }
  
  g <- ggscatter(fit, 
                x = "PC1", 
                y = "PC2",
                label = cpg_label,
                color = "groups",
                palette = "jco",
                font.label = 6,
                size = 2,
                ellipse = F,
                ellipse.type = "norm",
                repel = TRUE) + 
      xlab(pcalabs[1]) + 
      ylab(pcalabs[2]) +
      ggtitle(title) +
      theme_bw() + 
      theme(aspect.ratio=1,
            plot.title = element_text(hjust = 0.5))

  if(!is.null(output_file)){
    ggsave(output_file, g, width=5, height=4, units="in", dpi=300, bg="white")
  }  
  return(g)
}

cat("plot bvalue PCA of all CpG sites ...\n")
g = draw_pca_plot(as.matrix(cpg_bvalue_df), meta$group, paste0(project, ".CpG.all.PCA.png"))

cat("plot bvalue PCA of top 10000 variable CpG sites ...\n")
g = draw_pca_plot(as.matrix(top10000_cpg_bvalue_df), meta$group, paste0(project, ".CpG.top10000.PCA.png"))

#https://www.rdocumentation.org/packages/minfi/versions/1.18.4/topics/mdsPlot
#Use the mdsPlot function from minfi package to draw MDS plot by Euclidean distance.
draw_mds_plot <- function(cpg_bvalue_dist, groups=NULL, output_file=NULL){
  fit <- cmdscale(cpg_bvalue_dist) |> 
    as.data.frame() |>
    dplyr::rename("MDS_1" = 1, "MDS_2" = 2) |>
    tibble::rownames_to_column("Sample")

  hasGroup = !is.null(groups)
  if(hasGroup){
    fit$groups = groups
  }else{
    fit$groups = "All"
  }

  if(nrow(fit) > 20){
    cpg_label=NULL
  }else{
    cpg_label=fit$Sample
  }
  dms_plot <- ggscatter(fit, 
                        x = "MDS_1", 
                        y = "MDS_2",
                        label = cpg_label,
                        xlab = "MDS 1",
                        ylab = "MDS 2",
                        color = "groups",
                        palette = "jco",
                        font.label = 6,
                        size = 2,
                        ellipse = F,
                        ellipse.type = "norm",
                        repel = TRUE) + 
              theme_bw() + 
              theme(aspect.ratio=1)

  if(!is.null(output_file)){
    ggsave(output_file, dms_plot, width=5, height=4, units="in", dpi=300, bg="white")
  }  
  return(dms_plot)
}

draw_euclidean_mds_plot <- function(cpg_bvalue_df, groups, output_file=NULL){
  cpg_bvalue_dist <- dist(t(cpg_bvalue_df))
  return(draw_mds_plot(cpg_bvalue_dist, groups, output_file))
}

draw_corr_mds_plot <- function(cpg_bvalue_df, groups, output_file=NULL){
  cpg_bvalue_dist <- 1 - cor(cpg_bvalue_df, method = "pearson")
  return(draw_mds_plot(cpg_bvalue_dist, groups, output_file))
}

cat("plot bvalue MDS of all CpG sites ...\n")
g = draw_euclidean_mds_plot(as.matrix(cpg_bvalue_df), meta$group, paste0(project, ".CpG.euclidean_distance.all.MDS.png"))

cat("plot bvalue MDS of top 10000 variable CpG sites ...\n")
g = draw_euclidean_mds_plot(as.matrix(top10000_cpg_bvalue_df), meta$group, paste0(project, ".CpG.euclidean_distance.top10000.MDS.png"))

cat("plot bvalue MDS of all CpG sites ...\n")
g = draw_corr_mds_plot(as.matrix(cpg_bvalue_df), meta$group, paste0(project, ".CpG.pearson_corr.all.MDS.png"))

cat("plot bvalue MDS of top 10000 variable CpG sites ...\n")
g = draw_corr_mds_plot(as.matrix(top10000_cpg_bvalue_df), meta$group, paste0(project, ".CpG.pearson_corr.top10000.MDS.png"))

methykit_version <- as.character(packageVersion("methylKit"))
writeLines(paste0("Methykit,v", methykit_version), paste0(project, ".methylKit.version"))

