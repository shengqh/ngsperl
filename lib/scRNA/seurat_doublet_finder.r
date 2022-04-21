
source("scRNA_func.r")
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(DT)
library(data.table)
library(digest)
library(heatmap3)
library(cowplot)
library(scales)
library(stringr)
library(htmltools)
library(patchwork)
library(glmGamPoi)
library(DoubletFinder)

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
reduction<-myoptions$reduction
npcs<-as.numeric(myoptions$pca_dims)
doublet_rates<-as.numeric(myoptions$doublet_rates)

assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")
by_harmony<-reduction=="harmony"

prefix<-outFile

finalList=readRDS(parFile1)
obj<-finalList$obj

sweep.res <- paramSweep_v3(obj, PCs = 1:npcs, sct = by_sctransform)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn$pK=as.numeric(as.character(bcmvn$pK))

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

nExp <- round(ncol(obj) * doublet_rates)  # expect 4% doublets
obj <- doubletFinder_v3(obj, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:npcs)

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
options<-data.frame("name"=c("pN", "pK", "nExp", "npcs"),
"value"=c(0.25, pK_choose, nExp, npcs))
write.csv(options, paste0(outFile, ".options.csv"), row.names=F)

#ggplot(bcmvn, aes(pK, BCmetric)) + geom_point()


#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

