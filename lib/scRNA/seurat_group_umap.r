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


my_dimplot<-function(object, group.by, label, labels=NULL, title=NULL, scale_colors=NULL,random.seed=20220606,...){
  require(ggplot2)
  require(scales)
  g<-DimPlot(object, group.by=group.by, label=label, ...)
  if(!is.null(title)){
    g<-g+ggtitle(title)
  }
  if(!is.null(labels)){
    nclusters<-length(unlist(unique(object[[group.by]])))
    ccolors<-hue_pal()(nclusters)
    set.seed(random.seed)
    scolors<-sample(ccolors, size=nclusters)
    g<-g+scale_color_manual(values=scolors, labels = labels)
  }
  return(g)
}


for(label in c(TRUE, FALSE)){
  label_str=ifelse(label, ".label", ".nolabel")
  png(paste0("AG3669", ".all", label_str, ".umap.png"), width=3300, height=2000, res=300)
  g<-my_dimplot(object=obj, group.by="seurat_clusters", title = "AG3669", label=label, labels=levels(obj$seurat_cell_type))
  print(g)
  dev.off()
}
