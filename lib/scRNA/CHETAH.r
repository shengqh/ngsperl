
library(Seurat)
library(CHETAH)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

species=myoptions$species
prefix=myoptions$prefix

#load reference
load(parFile2)

finalList<-readRDS(parFile1)
obj<-finalList$obj

#keep umap only
obj@reductions$pca=NULL
obj@reductions$tsne=NULL

prefix=paste0(outFile, ".CHETAH")

obj.sce = as.SingleCellExperiment(obj)
obj.sce<-CHETAHclassifier(input=obj.sce, ref_cells = reference)

ct=data.frame(cell=c(1:length(colnames(obj.sce))), "seurat_cluster"=obj$seurat_clusters, "celltype_CHETAH"=obj.sce$celltype_CHETAH)
row.names(ct)<-colnames(obj)
write.csv(ct, paste0(prefix,".csv"))
saveRDS(obj.sce, file=paste0(prefix, ".rds"))

png(paste0(prefix, ".png"), width=6000, height=2500, res=300)
g<-PlotCHETAH(obj.sce)
print(g)
dev.off()

