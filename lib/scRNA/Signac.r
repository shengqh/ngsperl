library(SignacX)
library(Seurat)

finalList=readRDS(parFile1)
obj=finalList$obj

defaultAssay = DefaultAssay(obj)
if(defaultAssay == "integrated"){
  if("SCT" %in% names(obj@assays)){
    DefaultAssay(obj)="SCT"
  }else{
    DefaultAssay(obj)="RNA"
  }
  obj <- SCTransform(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:50, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:50, verbose = FALSE)  
}

labels <- SignacFast(obj)
celltypes_fast = GenerateLabels(labels, E = obj)

DefaultAssay(obj)=defaultAssay

obj <- AddMetaData(obj, metadata = celltypes_fast$CellStates, col.name = "CellStates")
png(paste0(outFile, ".CellStates.fast.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

saveRDS(celltypes_fast, file=paste0(outFile, ".CellStates.fast.rds"))

# labels <- Signac(obj)
# celltypes = GenerateLabels(labels, E = obj)

# obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "CellStates")
# png(paste0(outFile, ".CellStates.slow.png"), width=3300, height=3000, res=300)
# g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
# print(g1)
# dev.off()

# saveRDS(celltypes, file=paste0(outFile, ".CellStates.slow.rds"))
