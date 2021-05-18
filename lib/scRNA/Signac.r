library(SignacX)
library(Seurat)

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

obj <- AddMetaData(obj, metadata = celltypes_fast$CellStates, col.name = "CellStates")
png(paste0(outFile, ".signac.fast.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

saveRDS(celltypes_fast, file=paste0(outFile, ".signac.fast.rds"))

labels <- Signac(obj)
celltypes = GenerateLabels(labels, E = obj)

obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "CellStates")
png(paste0(outFile, ".signac.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

saveRDS(celltypes, file=paste0(outFile, ".signac.rds"))
