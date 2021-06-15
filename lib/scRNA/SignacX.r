rm(list=ls()) 
outFile='NoAAC_iSGS_airway_atlas'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/scratch/cqs/shengq2/alexander_gelbard_projects/20210415_NoAAC_iSGS_scRNA/seurat_harmony/result/NoAAC_iSGS_airway_atlas.final.rds'
parFile2=''
parFile3=''
tcell_only=0

setwd('C:/projects/scratch/cqs/shengq2/alexander_gelbard_projects/20210415_NoAAC_iSGS_scRNA/seurat_harmony_signac/result')

source("scRNA_func.r")
library(SignacX)
library(Seurat)
library(ggplot2)

finalList=readRDS(parFile1)
obj=finalList$obj

if(DefaultAssay(obj) == "integrated"){
  if(nrow(obj@assays$integrated@counts) == 0){
    if ("SCT" %in% names(obj@assays)){
      obj@assays$integrated@counts = obj@assays$SCT@counts
    }else{
      obj@assays$integrated@counts = obj@assays$RNA@counts
    }
  }
}

labels <- SignacFast(obj)
celltypes_fast = GenerateLabels(labels, E = obj)
saveRDS(celltypes_fast, file=paste0(outFile, ".SignacX.fast.rds"))

obj <- AddMetaData(obj, metadata = celltypes_fast$CellStates, col.name = "CellStates")
png(paste0(outFile, ".SignacX.fast.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

labels <- Signac(obj)
celltypes = GenerateLabels(labels, E = obj)
saveRDS(celltypes, file=paste0(outFile, ".SignacX.rds"))

obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "CellStates")
png(paste0(outFile, ".SignacX.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

dm=cbind(data.frame(obj@reductions$umap@cell.embeddings), data.frame(Cluster=obj$seurat_clusters, SignacX=obj$CellStates))
png(paste0(outFile, ".SignacX.cluster.png"), width=5000, height=4000, res=300)
g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=SignacX)) + geom_point() + facet_wrap(~Cluster)+theme_bw()
print(g1)
dev.off()

ct<-data.frame("Signac"=celltypes$CellStates, "Clusters"=obj$seurat_clusters)

ct_tbl<-table(ct$Signac,ct$Clusters)
write.csv(ct_tbl, paste0(outFile, ".SignacX.cluster.csv"))
