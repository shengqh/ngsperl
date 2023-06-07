rm(list=ls()) 
sample_name='DM_1'
outFile='DM_1'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/vickers_lab/projects/20230509_9061_scRNA_mouse_decontX_byTiger/decontX_doublet_finder/result/DM_1')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(DoubletFinder)

options(future.globals.maxSize= 10779361280)
random.seed=20200107
min.pct=0.5
logfc.threshold=0.6

myoptions<-read_file_map(parSampleFile2, do_unlist = FALSE)

by_sctransform<-is_one(myoptions$by_sctransform)
use_sctransform_v2<-is_one(myoptions$use_sctransform_v2)
npcs<-as.numeric(myoptions$pca_dims)
mc.cores<-as.numeric(myoptions$threads)

prefix<-outFile

if(file.exists(parFile1)){
  sample_file=parFile1
}else{
  file_map=read_file_map(parSampleFile1, do_unlist = FALSE)
  sample_file=file_map[[sample_name]]
}

lst=read_scrna_data(sample_file)
counts=lst$counts
obj = CreateSeuratObject(counts = counts, project = sample_name)
if(by_sctransform){
  obj=do_sctransform(obj, mc.cores=mc.cores, use_sctransform_v2=use_sctransform_v2)
  DefaultAssay(obj) <- "SCT"
}else{
  obj=do_normalization(obj)
  DefaultAssay(obj) <- "RNA"
}
obj = RunPCA(obj)
#obj = RunUMAP(obj)

sweep.res <- paramSweep_v3(obj, PCs = 1:npcs, sct = by_sctransform)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn$pK=as.numeric(as.character(bcmvn$pK))

pK=bcmvn$pK
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]
pN=0.25

#10x doublet rate
#https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
doublet_rates = seq(0.5, 1, 0.1) * ncol(obj) / 1000 * 0.008

params_df<-data.frame("pN"=pN, "pK"=pK_choose, "doublet_rate"=doublet_rates)
params_df$nExp = round(ncol(obj) * params_df$doublet_rate)
params_df$label=paste0("DF.classifications_", params_df$pN, "_", params_df$pK, "_", params_df$nExp)
write.csv(params_df, paste0(outFile, ".options.csv"), row.names=F)

doublet_rate=doublet_rates[1]
for (doublet_rate in doublet_rates){
  cat("doublet_rate=", doublet_rate, "\n")
  nExp <- round(ncol(obj) * doublet_rate)  # expect 4% doublets
  obj <- doubletFinder_v3(obj, pN =pN, pK = pK_choose, nExp = nExp, PCs = 1:npcs, sct=by_sctransform)
}

write.csv(obj@meta.data, paste0(outFile, ".meta.csv"))
saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
