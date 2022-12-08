rm(list=ls()) 
outFile='PL_7114_human'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='C:/projects/nobackup/h_turner_lab/shengq2/20220805_7114_scRNA_human/seurat_sct_harmony_multires_03_choose/result/PL_7114_human.final.rds'
parFile2=''
parFile3=''


setwd('C:/projects/nobackup/h_turner_lab/shengq2/20220805_7114_scRNA_human/seurat_sct_harmony_multires_04_group_umap/result')

### Parameter setting end ###

#library(scico)

source("scRNA_func.r")
if(!exists('obj')){
  obj<-read_object(parFile1)
}

groups<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
groups_map<-split(groups$V1, groups$V2)

#g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type")
#print(g)

heatlist=list()

label<-TRUE
for(label in c(TRUE, FALSE)){
  label_str=ifelse(label, ".label", ".nolabel")
  g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type")
  png(paste0(outFile, ".all", label_str, ".umap.png"), width=2500, height=2000, res=300)
  print(g)
  dev.off()

  gname = names(groups_map)[1]
  for(gname in names(groups_map)){
    samples = groups_map[[gname]]
    cells<-colnames(obj)[obj$orig.ident %in% samples]
    subobj<-subset(obj, cells=cells)
    g<-get_dim_plot(subobj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = gname, legend.title = "Cell type")
    png(paste0(outFile, ".", gname, label_str, ".umap.png"), width=2500, height=2000, res=300)
    print(g)
    dev.off()
    
    if(label){
      coord=FetchData(obj, c("UMAP_1", "UMAP_2", "orig.ident"))
      in_group<-coord[colnames(subobj),]
      out_group<-coord[setdiff(rownames(coord), colnames(subobj)),]
      
      g<-ggplot(out_group, aes(UMAP_1, UMAP_2)) + geom_point(color="grey") + 
        geom_hex(aes(UMAP_1, UMAP_2), data=in_group, bins=70) +
        scale_fill_gradient(low="white", high="red") + ggtitle(gname) + theme_bw3()
      
      heatlist[[gname]] = g
    }
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

library(patchwork)

nf<-length(heatlist)
ncol=ceiling(sqrt(nf))
nrow=ceiling(nf/ncol)

gg<-patchwork::wrap_plots(heatlist, nrow =nrow, ncol=ncol)
png(paste0(outFile, ".heat.png"),width=ncol * 1500, height=nrow*1300, res=300)
print(gg)
dev.off()

save_highlight_cell_plot(paste0(outFile, ".cell.png"), obj, group.by = "seurat_cell_type")

