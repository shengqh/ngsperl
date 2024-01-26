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

meta <- read.table(parSampleFile3, sep = "\t", header = F)
colnames(meta) <- c("sample", "group")
var = "group"

if(!is.na(control_group)){
  stopifnot(any(meta$group == control_group))
  treatment<-ifelse(meta$group == control_group, 0, 1)
}else{
  treatment<-NA
}

rds_file=paste0(project, ".filtered.cpg.meth.rds")
#if(!file.exists(rds_file)){
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

  saveRDS(filtered.cpg.meth, rds_file)
#}else{
#  filtered.cpg.meth<-readRDS(rds_file)
#}

cpg_bvalue_df <- data.frame(filtered.cpg.meth) %>%
  mutate(id = paste(chr, start, end, sep = "_"))

sample.ids=filtered.cpg.meth@sample.ids
for(i in 1:length(sample.ids)){
  id <- sample.ids[i]
  ncs <- paste0("numCs", i)
  cov <- paste0("coverage", i)
  cpg_bvalue_df[, id] <- cpg_bvalue_df[, ncs] / cpg_bvalue_df[, cov]
}
cpg_bvalue_df <- cpg_bvalue_df[,c("id", file.id)] %>%
  column_to_rownames("id")
saveRDS(cpg_bvalue_df, paste0(project, ".cpg_bvalue_df.rds"))

cpg_bvalue_cor <- cor(cpg_bvalue_df, method = "pearson")
saveRDS(cpg_bvalue_cor, paste0(project, ".cpg_bvalue_cor.rds"))

cpg_bvalue_mds <- (1 - cpg_bvalue_cor) %>%
  cmdscale() %>%
  data.frame()
colnames(cpg_bvalue_mds) <- c("Dim.1", "Dim.2")

stopifnot(rownames(cpg_bvalue_mds) == colnames(cpg_bvalue_df))

cpg_bvalue_mds<-cpg_bvalue_mds[meta$sample,]
stopifnot(rownames(cpg_bvalue_mds) == meta$sample)

cpg_bvalue_mds[,var] <- meta[,var]
saveRDS(cpg_bvalue_mds, paste0(project, ".cpg_bvalue_mds.rds"))

# Plot MDS
if(nrow(cpg_bvalue_mds) > 20){
  cpg_label=NULL
}else{
  cpg_label=rownames(cpg_bvalue_mds)
}
dms_plot <- ggscatter(cpg_bvalue_mds, x = "Dim.1", y = "Dim.2",
                      label = cpg_label,
                      xlab = "MDS 1",
                      ylab = "MDS 2",
                      color = var,
                      palette = "jco",
                      font.label = 2,
                      size = 1,
                      ellipse = F,
                      ellipse.type = "norm",
                      repel = TRUE) + theme_bw()

ggsave(paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.png"), dms_plot, width=4, height=3, units="in", dpi=300, bg="white")
ggsave(paste0(project, "_methyl_CpG_bvalue_corr_MDS_plot.pdf"), dms_plot, width=4, height=3, units="in", dpi=300, bg="white")

