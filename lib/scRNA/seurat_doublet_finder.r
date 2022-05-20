
library(Seurat)
library(DoubletFinder)

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
npcs<-as.numeric(myoptions$pca_dims)

prefix<-outFile

if(!exists("obj")){
  obj=read_object(parFile1, parFile2)
}

sweep.res <- paramSweep_v3(obj, PCs = 1:npcs, sct = by_sctransform)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn$pK=as.numeric(as.character(bcmvn$pK))

pK=bcmvn$pK
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pN=0.25

params_df<-data.frame("pN"=pN, "pK"=pK_choose, "doublet_rate"=c(1:10)*0.01)
params_df$nExp = round(ncol(obj) * params_df$doublet_rate)
params_df$label=paste0("DF.classifications_", params_df$pN, "_", params_df$pK, "_", params_df$nExp)
write.csv(params_df, paste0(outFile, ".options.csv"), row.names=F)

doublet_rate=0.01
for (doublet_rate in c(1:10) * 0.01){
  cat("doublet_rate=", doublet_rate, "\n")
  nExp <- round(ncol(obj) * doublet_rate)  # expect 4% doublets
  obj <- doubletFinder_v3(obj, pN =pN, pK = pK_choose, nExp = nExp, PCs = 1:npcs, sct=by_sctransform)
}

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
  
