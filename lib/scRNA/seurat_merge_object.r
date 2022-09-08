rm(list=ls()) 
outFile='PH_combine'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('C:/projects/scratch/cqs/shengq2/paula_hurley_projects/20220824_scRNA_7467_benign_hg38/seurat_merge_object/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(digest)
library(patchwork)
library(sparseMatrixStats)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile2, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

sample_pattern= myoptions$sample_pattern

Mtpattern= myoptions$Mtpattern
rRNApattern=myoptions$rRNApattern
hemoglobinPattern=myoptions$hemoglobinPattern


#read raw count dat
filelist1<-read.table(parSampleFile1, header=F, stringsAsFactors = F)
rawobjs = list()
fidx=1
fileMap<-split(filelist1$V1, filelist1$V2)
fileTitle<-names(fileMap)[1]

for(fileTitle in names(fileMap)) {
  fileName  = fileMap[[fileTitle]]
  cat(fileTitle, "\n")
  sobj<-readRDS(fileName)
  if(is.list(sobj)){
    sobj<-sobj$rawobj
  }
  
  if(sample_pattern != ""){
    cells = colnames(sobj)[grepl(sample_pattern, sobj$orig.ident)]
    sobj<-subset(sobj, cells=cells)
  }

  rawobjs[[fileTitle]] = sobj
  rm(sobj)
}

if(length(rawobjs) == 1){
  rawobj <- rawobjs[[1]]
}else{
  rawobj <- merge(rawobjs[[1]], y = unlist(rawobjs[2:length(rawobjs)]), project = "integrated")
}
rm(rawobjs)

writeLines(rownames(rawobj), paste0(outFile, ".genes.txt"))

rawobj<-PercentageFeatureSet(object=rawobj, pattern=Mtpattern, col.name="percent.mt")
rawobj<-PercentageFeatureSet(object=rawobj, pattern=rRNApattern, col.name = "percent.ribo")
rawobj<-PercentageFeatureSet(object=rawobj, pattern=hemoglobinPattern, col.name="percent.hb")    

saveRDS(rawobj, paste0(outFile, ".rawobj.rds"))

png(paste0(outFile, ".top20.png"), width=3000, height=2000, res=300)
par(mar = c(4, 8, 2, 1))
C <- rawobj@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
mc<-rowMedians(C)
most_expressed <- order(mc, decreasing = T)[20:1]
tm<-as.matrix(Matrix::t(C[most_expressed,]))
boxplot(tm, cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()

draw_feature_qc(outFile, rawobj, "orig.ident")

if(any(rawobj$orig.ident != rawobj$sample)){
  draw_feature_qc(paste0(outFile, ".sample"), rawobj, "sample")
}

rRNA.genes <- grep(pattern = rRNApattern,  rownames(rawobj), value = TRUE)
rawobj<-rawobj[!(rownames(rawobj) %in% rRNA.genes),]

rawobj<-PercentageFeatureSet(object=rawobj, pattern=Mtpattern, col.name="percent.mt")
rawobj<-PercentageFeatureSet(object=rawobj, pattern=rRNApattern, col.name = "percent.ribo")
rawobj<-PercentageFeatureSet(object=rawobj, pattern=hemoglobinPattern, col.name="percent.hb")    

draw_feature_qc(paste0(outFile, ".no_ribo"), rawobj, "orig.ident")

if(any(rawobj$orig.ident != rawobj$sample)){
  draw_feature_qc(paste0(outFile, ".no_ribo", ".sample"), rawobj, "sample")
}

