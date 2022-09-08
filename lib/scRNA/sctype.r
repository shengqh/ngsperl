library(Seurat)
library(ggplot2)
library("dplyr")
library("HGNChelper")

options(future.globals.maxSize= 10779361280)
random.seed=20200107

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)
tcell_only = ifelse(myoptions$tcell_only == "1", TRUE, FALSE)
pca_dims=1:as.numeric(myoptions$pca_dims)
resolution=as.numeric(myoptions$resolution)
reduction=myoptions$reduction
tissue=myoptions$tissue

if(!exists("obj")){
  obj=read_object(parFile1, parFile2)
}

if(tcell_only){
  cells_tbl = obj@meta.data
  cells=rownames(cells_tbl)[cells_tbl[,myoptions$celltype_layer] == "T cells"]
  obj=subset(obj, cells=cells)
}

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", tissue) # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

obj=ScaleData(obj, features = rownames(obj))

es.max = sctype_score(scRNAseqData = obj[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(obj@meta.data[obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
#sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
sct=split(sctype_scores$type, sctype_scores$cluster)
obj$sctype=unlist(sct[as.character(obj$seurat_clusters)])

saveRDS(obj@meta.data, file=)

saveRDS(celltypes, file=paste0(outFile, ".SignacX.rds"))
