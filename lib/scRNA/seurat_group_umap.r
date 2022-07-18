rm(list=ls()) 
outFile='AK6383'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_03_choose_res/result/AK6383.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/nobackup/kirabo_lab/shengq2/20220506_6383_scRNA_human/seurat_merge_04_group_umap/result')

### Parameter setting end ###

source("scRNA_func.r")
obj<-read_object(parFile1)

groups<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
groups_map<-unlist(split(groups$V2, groups$V1))

obj$group <- groups_map[obj$orig.ident]

g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type")
print(g)

for(label in c(TRUE, FALSE)){
  label_str=ifelse(label, ".label", ".nolabel")
  g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type")
  png(paste0(outFile, ".all", label_str, ".umap.png"), width=2500, height=2000, res=300)
  print(g)
  dev.off()
  
  ngroups=length(unique(obj$group))
  g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type", split.by="group")
  g<-add_x_y_axis(g)
  png(paste0(outFile, ".group", label_str, ".umap.png"), width=2000 * ngroups + 500, height=2000, res=300)
  print(g)
  dev.off()
}

save_highlight_cell_plot(paste0(outFile, ".cell.png"), obj, group.by = "seurat_cell_type")

