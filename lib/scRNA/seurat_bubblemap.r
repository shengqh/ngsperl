source("scRNA_func.r")

library("Seurat")
library("readxl")
library(ggplot2)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")

finalList<-readRDS(parFile1)
obj=finalList$obj

genes <- read_xlsx(parFile4, sheet = 1)
for(idx in c(2:nrow(genes))){
  if(is.na(genes[idx,"Cell Type"])){
    genes[idx,"Cell Type"]=genes[idx-1,"Cell Type"]
  }
}
genes$`Cell Type`=factor(genes$`Cell Type`, levels=unique(genes$`Cell Type`))

gene_names=genes$`Marker Gene`
gene_names[gene_names=="PECAM"] = "PECAM1"
gene_names[gene_names=="HGD1B"] = "HGD"
gene_names[gene_names=="EpCAM"] = "EPCAM"
gene_names=factor(gene_names, levels=gene_names)
genes$`Marker Gene`<-gene_names

gene_groups=split(genes$`Marker Gene`, genes$`Cell Type`)

cell_type<-read.csv(parFile3)

sheets=excel_sheets(parFile4)
if(length(sheets) > 1){
  clusters<-read_xlsx(parFile4, sheet = 2)
  cluster_ids<-clusters$`Order of Clusters`
  cell_type$seurat_clusters<-factor(cell_type$seurat_clusters, levels=cluster_ids)
  cell_type<-cell_type[order(cell_type$seurat_clusters),]
  cell_type$seurat_renamed_cellactivity_clusters=factor(cell_type$seurat_renamed_cellactivity_clusters, levels=(cell_type$seurat_renamed_cellactivity_clusters))
  rownames(cell_type)=cell_type$seurat_clusters
  obj[["seurat_renamed_cellactivity_clusters"]]=cell_type[as.character(obj$seurat_clusters),"seurat_renamed_cellactivity_clusters"]
  group.by="seurat_renamed_cellactivity_clusters"
}else{
  ct<-cell_type[!duplicated(cell_type$cell_type),]
  cell_type$cell_type<-factor(cell_type$cell_type, levels=ct$cell_type)
  cell_type<-cell_type[order(cell_type$cell_type, cell_type$seurat_clusters),]
  cell_type$seurat_celltype_clusters=paste0(cell_type$seurat_clusters, " : ", cell_type$cell_type)
  cell_type$seurat_celltype_clusters=factor(cell_type$seurat_celltype_clusters, levels=cell_type$seurat_celltype_clusters)
  rownames(cell_type)=cell_type$seurat_clusters
  obj[["seurat_celltype_clusters"]]=cell_type[as.character(obj$seurat_clusters),"seurat_celltype_clusters"]
  group.by="seurat_celltype_clusters"
}

png(paste0(outFile, ".bubblemap.png"), width=5500, height=3000, res=300)
g=DotPlot(obj, features=gene_groups, assay=assay, group.by=group.by) + 
  xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                              strip.background = element_blank(),
                              strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
print(g)
dev.off()
