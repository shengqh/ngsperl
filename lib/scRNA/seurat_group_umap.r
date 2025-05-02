rm(list=ls()) 
outFile='Q51804_liver'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1='/nobackup/h_cqs/zhaos/Patterson/202504_snRNA_Q51804_liver/seurat_sct2_fastmnn_dr0.5_3_choose/result/Q51804_liver.final.rds'
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/zhaos/Patterson/202504_snRNA_Q51804_liver/seurat_sct2_fastmnn_dr0.5_3_choose_group_umap/result')

### Parameter setting end ###

source("scRNA_func.r")
#library(scico)

obj<-read_object(parFile1)

groups<-read.table(parSampleFile1, sep="\t", stringsAsFactors = F)
groups_map<-split(groups$V1, groups$V2)

group_names=names(groups_map)
ngroup=length(group_names)

#g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type")
#print(g)

heatlist=list()

save_plot=function(filepath, g){
  nclusters=length(unique(g$seurat_clusters))
  if(nclusters > 15){
    g=g+guides(fill = guide_legend(ncol = 2))
    width=14
  }else{
    width=11
  }
  height=6
  ggsave(filepath, g, width=width, height=height, units="in", dpi=300, bg="white")
}

label<-TRUE
for(label in c(TRUE, FALSE)){
  label_str=ifelse(label, ".label", ".nolabel")
  g<-get_dim_plot(obj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type") +
    theme(plot.title=element_blank()) +
    xlab("UMAP_1") + ylab("UMAP_2")
  save_plot(paste0(outFile, ".All", label_str, ".umap.png"), g)

  gname = group_names[2]
  for(gname in group_names){
    if(gname == "All"){
      next
    }
    cat("Processing group: ", gname, "\n")

    samples = groups_map[[gname]]
    cells<-colnames(obj)[obj$orig.ident %in% samples]
    subobj<-subset(obj, cells=cells)

    g<-get_dim_plot(subobj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = gname, legend.title = "Cell type") +
      theme(plot.title=element_blank()) +
      xlab("UMAP_1") + ylab("UMAP_2")
    save_plot(paste0(outFile, ".", gname, label_str, ".umap.png"), g)
    
    if(label){
      g=DimPlot(obj, group.by="orig.ident")
      coord=g$data |> dplyr::rename(UMAP_1=1, UMAP_2=2)
      in_group<-coord[colnames(subobj),]
      out_group<-coord[setdiff(rownames(coord), colnames(subobj)),]
      
      g<-ggplot(out_group, aes(UMAP_1, UMAP_2)) + geom_point(color="grey") + 
        geom_hex(aes(UMAP_1, UMAP_2), data=in_group, bins=70) +
        scale_fill_gradient(low="white", high="red") + ggtitle(gname) + theme_bw3()
      
      heatlist[[gname]] = g
    }
  }

  if(ngroup > 1){
    if(all(table(groups$V1) == 1)){ 
      cat("Processing multiple group figure ...\n")
      #each sample belongs to one group
      cells<-colnames(obj)[obj$orig.ident %in% groups$V1]
      subobj<-subset(obj, cells=cells)
      cur_groups_map<-unlist(split(groups$V2, groups$V1))
      subobj@meta.data$group <- cur_groups_map[subobj$orig.ident]

      ngroups=length(unique(subobj$group))
      g<-get_dim_plot(subobj, group.by="seurat_clusters", label.by = "seurat_cell_type", label=label, title = "", legend.title = "Cell type", split.by="group")
      g<-add_x_y_axis(g)
      ggsave(paste0(outFile, ".group", label_str, ".umap.png"), g,  width=2000 * ngroups + 500, height=2000, units="px", dpi=300, bg="white")
    }
  }
}

if(ngroup > 1){
  library(patchwork)

  nf<-length(heatlist)
  ncol=ceiling(sqrt(nf))
  nrow=ceiling(nf/ncol)

  gg<-patchwork::wrap_plots(heatlist, nrow =nrow, ncol=ncol)
  ggsave(paste0(outFile, ".heat.png"), gg, width=ncol * 1500, height=nrow*1300, units="px", dpi=300, bg="white")

  save_highlight_cell_plot(paste0(outFile, ".cell.png"), obj, group.by = "seurat_cell_type")
}

