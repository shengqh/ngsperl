library(SignacX)
library(Seurat)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
tcell_only = ifelse(myoptions$tcell_only == "1", TRUE, FALSE)
pca_dims=1:as.numeric(myoptions$pca_dims)
resolution=as.numeric(myoptions$resolution)
reduction=myoptions$reduction
by_sctransform<-ifelse(myoptions$by_sctransform == "0", FALSE, TRUE)
assay=ifelse(by_sctransform, "SCT", "RNA")


if(!exists("obj")){
  obj=read_object(parFile1, parFile2)
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
  cells_tbl = obj@meta.data
  cells=rownames(cells_tbl)[cells_tbl[,myoptions$celltype_layer] == "T cells"]
  obj=subset(obj, cells=cells)
}

obj<-FindNeighbors(object = obj, reduction=reduction, dims=pca_dims, verbose=FALSE)

labels <- Signac(obj)

celltypes = GenerateLabels(labels, E = obj)
saveRDS(celltypes, file=paste0(outFile, ".SignacX.rds"))

obj <- AddMetaData(obj, metadata = celltypes$CellStates, col.name = "signacx_CellStates")

bubblemap_file=myoptions$bubblemap_file
has_bubblemap <- !is.null(bubblemap_file) && file.exists(bubblemap_file)

g1=DimPlot(obj, group.by = "signacx_CellStates", reduction="umap", label=T)
g2=DimPlot(obj, group.by = "seurat_clusters", reduction="umap", label=T)

if(has_bubblemap){
  g3<-get_bubble_plot(obj, NA, "signacx_CellStates", bubblemap_file, assay="RNA")
  layout <- "AB
CC"
  g<-g1+g2+g3+plot_layout(design=layout)
  height=6000
}else{
  g<-g1+g2+plot_layout(ncol=2)
  height=3000
}

png(paste0(outFile, ".SignacX.no_recluster.png"), width=6300, height=height, res=300)
print(g)
dev.off()

dm=cbind(data.frame(obj@reductions$umap@cell.embeddings), data.frame(Cluster=obj$seurat_clusters, SignacX=obj$signacx_CellStates))
dm$Cluster=factor(dm$Cluster)

png(paste0(outFile, ".SignacX.no_recluster.c1.png"), width=5000, height=4000, res=300)
g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=SignacX)) + geom_point() + facet_wrap(~Cluster)+theme_bw()
print(g1)
dev.off()

png(paste0(outFile, ".SignacX.no_recluster.c2.png"), width=5000, height=4000, res=300)
g1=ggplot(dm, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) + geom_point() + facet_wrap(~SignacX)+theme_bw()
print(g1)
dev.off()

ct<-data.frame("SignacX"=obj$signacx_CellStates, "Clusters"=obj$seurat_clusters)

ct_tbl<-table(ct$Signac,ct$Clusters)
write.csv(ct_tbl, paste0(outFile, ".SignacX.cluster.csv"))

saveRDS(obj@meta.data, paste0(outFile, ".meta.rds"))

