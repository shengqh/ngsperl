rm(list=ls()) 
outFile='AG_integrated'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/data/h_gelbard_lab/projects/20220907_8566_project/seurat_sct_harmony_multires_03_choose/result/AG_integrated.final.rds'
parFile2=''
parFile3=''


setwd('/data/h_gelbard_lab/projects/20220907_8566_project/seurat_sct_harmony_multires_04_group_umap/result')

### Parameter setting end ###

source("scRNA_func.r")
if(!exists('obj')){
  obj<-read_object(parFile1)
}

g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type")
print(g)

label<-TRUE
for(label in c(TRUE, FALSE)){
  label_str=ifelse(label, ".label", ".nolabel")
  g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type")
  png(paste0(outFile, ".all", label_str, ".umap.png"), width=2500, height=2000, res=300)
  print(g)
  dev.off()

  groups<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
  groups_map<-split(groups$V1, groups$V2)

  for(gname in names(groups_map)){
    samples = groups_map[[gname]]
    cells<-colnames(obj)[obj$orig.ident %in% samples]
    subobj<-subset(obj, cells=cells)
    g<-get_dim_plot(subobj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = gname, legend.title = "Cell type")
    png(paste0(outFile, ".", gname, label_str, ".umap.png"), width=2500, height=2000, res=300)
    print(g)
    dev.off()
  }

  if(all(table(groups$V1) == 1)){ 
    #each sample belongs to one group
    cells<-colnames(obj)[obj$orig.ident %in% groups$V1]
    subobj<-subset(obj, cells=cells)
    groups_map<-unlist(split(groups$V2, groups$V1))
    subobj$group <- groups_map[subobj$orig.ident]

    ngroups=length(unique(subobj$group))
    g<-get_dim_plot(subobj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type", split.by="group")
    g<-add_x_y_axis(g)
    png(paste0(outFile, ".group", label_str, ".umap.png"), width=2000 * ngroups + 500, height=2000, res=300)
    print(g)
    dev.off()
  }
}

save_highlight_cell_plot(paste0(outFile, ".cell.png"), obj, group.by = "seurat_cell_type")

