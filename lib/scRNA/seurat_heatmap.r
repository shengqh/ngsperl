source("scRNA_func.r")

library("Seurat")
library("readxl")
library(ggplot2)

options_table<-read.table(parSampleFile1, sep="\t", header=F, stringsAsFactors = F)
myoptions<-split(options_table$V1, options_table$V2)

assay=ifelse(myoptions$by_sctransform == "0", "RNA", "SCT")

finalList<-readRDS(parFile1)
obj=finalList$obj

assaydata=GetAssayData(obj, assay=assay)
allgenes=rownames(assaydata)
rm(assaydata)

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
gene_names[gene_names=="CD25"] = "IL2RA"
gene_names[gene_names=="ACTAA2"] = "ACTA2"
gene_names[gene_names=="MTND6"] = "MT-ND6"
gene_names[gene_names=="FOXJ!"] = "FOXJ1"

genes$`Marker Gene`<-gene_names

miss_genes=genes$`Marker Gene`[!(genes$`Marker Gene` %in% allgenes)]
writeLines(miss_genes, con="miss_gene.csv")

genes<-genes[genes$`Marker Gene` %in% allgenes,]

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
  ct_levels<-c("B cells", "Plasma cells", "NK cells", "T cells", "Macrophages", "Mast cells", "Endothelial cells", "Fibroblasts", "Epithelial cells", "Basal cells", "Olfactory epithelial cells", "Ciliated cells")
  ct<-cell_type[!duplicated(cell_type$cell_type),]
  missed = ct$cell_type[!(ct$cell_type %in% ct_levels)]
  if(length(missed) > 0){
    ct_levels = c(ct_levels, missed)
  }
  ct_levels = ct_levels[ct_levels %in% ct$cell_type]
  cell_type$cell_type<-factor(cell_type$cell_type, levels=ct_levels)
  cell_type<-cell_type[order(cell_type$cell_type, cell_type$seurat_clusters),]
  cell_type$seurat_celltype_clusters=paste0(cell_type$seurat_clusters, " : ", cell_type$cell_type)
  cell_type$seurat_celltype_clusters=factor(cell_type$seurat_celltype_clusters, levels=cell_type$seurat_celltype_clusters)
  rownames(cell_type)=cell_type$seurat_clusters
  obj[["seurat_celltype_clusters"]]=cell_type[as.character(obj$seurat_clusters),"seurat_celltype_clusters"]
  group.by="seurat_celltype_clusters"
}

genes=unique(unlist(gene_groups))
g<-DotPlot(obj, features=genes, assay="RNA",group.by=group.by)
gdata<-g$data

data.plot<-NULL
gn=names(gene_groups)[1]
for(gn in names(gene_groups)){
  gs=gene_groups[[gn]]
  gdd<-gdata[gdata$features.plot %in% gs,]
  if(nrow(gdd)== 0){
    stop(gn)
  }
  gdd$feature.groups=gn
  data.plot<-rbind(data.plot, gdd)
}

data.plot$feature.groups=factor(data.plot$feature.groups, levels=names(gene_groups))

color.by <- "avg.exp.scaled"
scale.func <- scale_radius
scale.min = NA
scale.max = NA
dot.scale = 6
cols = c("lightgrey", "blue")

library(cowplot)
plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
  geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
  scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
  labs(x = "Features", y = "Identity") +
  theme_cowplot() + 
  facet_grid(facets = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") + 
  theme(panel.spacing = unit(x = 1,units = "lines"), strip.background = element_blank()) + 
  scale_color_gradient(low = cols[1], high = cols[2])


png(paste0(outFile, ".bubblemap.png"), width=5500, height=3000, res=300)
g=plot + 
  xlab("") + ylab("") + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
                              strip.background = element_blank(),
                              strip.text.x = element_text(angle=90, hjust=0, vjust=0.5))
print(g)
dev.off()
