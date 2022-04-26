library(SignacX)
library(Seurat)
library(ggplot2)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
tcell_only = ifelse(myoptions$tcell_only == "1", TRUE, FALSE)
pca_dims=1:as.numeric(myoptions$pca_dims)
resolution=as.numeric(myoptions$resolution)
reduction=myoptions$reduction

if(!exists('obj')){
  obj=readRDS(parFile1)
  if(is.list(obj)){
    obj=obj$obj
  }
}

if(DefaultAssay(obj) == "integrated"){
  if(nrow(obj@assays$integrated@counts) == 0){
    if ("SCT" %in% names(obj@assays)){
      obj@assays$integrated@counts = obj@assays$SCT@counts
    }else{
      obj@assays$integrated@counts = obj@assays$RNA@counts
    }
  }
}

if(tcell_only){
  cells_tbl = read.csv(parFile2, row.names=1)
  cells=rownames(cells_tbl)[cells_tbl[,myoptions$celltype_layer] == "T cells"]
  obj=subset(obj, cells=cells)
  obj<-FindNeighbors(object = obj, reduction=reduction, dims=pca_dims, verbose=FALSE)
}

labels <- Signac(obj)

celltypes = GenerateLabels(labels, E = obj)
saveRDS(celltypes, file=paste0(outFile, ".SignacX.rds"))

obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "CellStates")
png(paste0(outFile, ".SignacX.no_recluster.png"), width=3300, height=3000, res=300)
g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
print(g1)
dev.off()

dm=cbind(data.frame(obj@reductions$umap@cell.embeddings), data.frame(Cluster=obj$seurat_clusters, SignacX=obj$CellStates))
png(paste0(outFile, ".SignacX.no_recluster.c1.png"), width=5000, height=4000, res=300)
g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=SignacX)) + geom_point() + facet_wrap(~Cluster)+theme_bw()
print(g1)
dev.off()

png(paste0(outFile, ".SignacX.no_recluster.c2.png"), width=5000, height=4000, res=300)
g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) + geom_point() + facet_wrap(~SignacX)+theme_bw()
print(g1)
dev.off()

ct<-data.frame("Signac"=celltypes$CellStates, "Clusters"=obj$seurat_clusters)

ct_tbl<-table(ct$Signac,ct$Clusters)
write.csv(ct_tbl, paste0(outFile, ".SignacX.cluster.csv"))

if(tcell_only){
  obj<-RunUMAP(object = obj, reduction=reduction, dims=pca_dims, verbose = FALSE)
  obj<-FindClusters(object=obj, verbose=FALSE, random.seed=random.seed, resolution=resolution)
  
  png(paste0(outFile, ".SignacX.recluster.png"), width=3300, height=3000, res=300)
  g1=DimPlot(obj, group.by = "CellStates", reduction="umap", label=T)
  print(g1)
  dev.off()
  
  dm=cbind(data.frame(obj@reductions$umap@cell.embeddings), data.frame(Cluster=obj$seurat_clusters, SignacX=obj$CellStates))
  png(paste0(outFile, ".SignacX.recluster.c1.png"), width=5000, height=4000, res=300)
  g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=SignacX)) + geom_point() + facet_wrap(~Cluster)+theme_bw()
  print(g1)
  dev.off()
  
  png(paste0(outFile, ".SignacX.recluster.c2.png"), width=5000, height=4000, res=300)
  g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) + geom_point() + facet_wrap(~SignacX)+theme_bw()
  print(g1)
  dev.off()
}

