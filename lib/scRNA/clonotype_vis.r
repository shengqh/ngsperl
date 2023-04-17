rm(list=ls()) 
outFile='AG3669'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20220508_scRNA_3669/seurat_merge_03_choose_res/result/AG3669.final.rds'
parFile2=''
parFile3=''
outputDirectory='.'


setwd('/data/h_gelbard_lab/projects/20220508_scRNA_3669/clonotype_03_convert_vis/result')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(ggplot2)
library(ggpubr)

obj<-read_object(parFile1)

if(parFile2 != ""){
  cell_df<-read_cell_cluster_file(parFile2)
}else{
  cell_df=obj@meta.data
}

if("seurat_cell_type" %in% colnames(cell_df)){
  cell_df$seurat_cellactivity_clusters = cell_df$seurat_cell_type
}

file_map = read.table("fileList1.txt", sep="\t")
clonotypes=NULL
for(parFile3 in file_map$V1){
  cur_clonotypes<-read.csv(parFile3)
  cur_clonotypes<-cur_clonotypes[order(cur_clonotypes$frequency, cur_clonotypes$clonotype_id, decreasing = T),]
  cur_clonotypes<-cur_clonotypes[c(1:min(10, nrow(cur_clonotypes))),,drop=F]
  clonotypes = rbind(clonotypes, cur_clonotypes)
}

clonocells<-unlist(strsplit(clonotypes$cells, split = ";"))
cell_df_clono<-subset(cell_df, rownames(cell_df) %in% clonocells)
cell_df_clono<-cell_df_clono[order(cell_df_clono$seurat_clusters),]
clono_clusters<-unique(cell_df_clono$seurat_clusters)
display_clusters<-as.character(unique(cell_df_clono$seurat_cellactivity_clusters))

valid_cell_df<-subset(cell_df, cell_df$seurat_clusters %in% clono_clusters)

valid_obj=subset(obj, cells=rownames(valid_cell_df))
rm(obj)

valid_obj$final_seurat_clusters=valid_cell_df$seurat_cellactivity_clusters

gcell<-MyDimPlot(valid_obj, group.by="final_seurat_clusters", reduction = "umap", label=T) + theme(legend.position = "none") + ggtitle("")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

idx=1
for(idx in c(1:nrow(clonotypes))) {
  curname=clonotypes$clonotype_id[idx]
  cdr3s_aa=clonotypes$cdr3s_aa[idx]
  cdr3s_cells=strsplit(clonotypes$cells[idx], split = ";")
  names(cdr3s_cells)<-cdr3s_aa
  cg<-MyDimPlot(valid_obj, reduction = "umap", cells.highlight = cdr3s_cells, label=F, repel=T) + theme(legend.position = "none") + 
    scale_color_manual(values = c("lightgrey", "red")) + ggtitle(cdr3s_aa)
  
  cobj=subset(valid_obj, cells=unlist(cdr3s_cells))
  sc<-data.frame(cell=colnames(cobj), sample=cobj$orig.ident)
  samplelist<-tapply(sc$cell,sc$sample,list)
  sg<-MyDimPlot(valid_obj, reduction = "umap", cells.highlight = samplelist, cols.highlight=gg_color_hue(length(samplelist)) , label=F) + theme(legend.position = "top")
  
  pdf(file=paste0(curname, ".umap.pdf"), width=21, height=7)
  p<-ggarrange(gcell, cg, sg, ncol=3)
  print(p)
  dev.off()
}
