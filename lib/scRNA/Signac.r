library(SignacX)
library(Seurat)

finalList=readRDS(parFile1)
obj=finalList$obj

labels <- SignacFast(obj)
celltypes_fast = GenerateLabels(labels, E = obj)

obj <- AddMetaData(obj, metadata = celltypes_fast$CellStates, col.name = "CellStates")
png(paste0(outFile, ".CellStates.fast.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

saveRDS(celltypes_fast, file=paste0(outFile, ".CellStates.fast.rds"))

labels <- Signac(obj)
celltypes = GenerateLabels(labels, E = obj)

obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "CellStates")
png(paste0(outFile, "CellStates.slow.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

saveRDS(celltypes, file=paste0(outFile, ".CellStates.slow.rds"))
