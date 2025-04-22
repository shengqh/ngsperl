rm(list=ls()) 
sample_name='WHY_01'
outFile='WHY_01'
parSampleFile1='fileList1.txt'
parSampleFile2=''
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/temp/20250414_T01_prepare_bin_counts/spatial_multi_umap/result/WHY_01')

### Parameter setting end ###

source("scRNA_func.r")
set.seed(20250414)

library(Seurat)

plan("multicore", workers = 16)
options(future.globals.maxSize = 10000 * 1024^2)
options(future.seed = TRUE)

file_map=read_file_map(parSampleFile1, do_unlist = FALSE)
data_path = file_map[[sample_name]]

obj=read_object(data_path, sample_name=sample_name)

df=NULL

n.neighbors=10
for(n.neighbors in c(10, 20, 30)) {
  min.dist=0.1
  for(min.dist in c(0.1, 0.2, 0.3)) {
    reduction.key=paste0("UMAP_N",n.neighbors,"_D",min.dist)
    # obj = RunUMAP(obj, 
    #   dims = 1:30,
    #   reduction = "pca",
    #   n.neighbors = n.neighbors,
    #   min.dist = min.dist)
    # g = get_dim_plot_labelby(obj, reduction = "umap", label.by = "cell_type") + 
    #   ggtitle(reduction.key)
    umap_png = paste0(sample_name, ".", reduction.key, ".png")
    # ggsave(umap_png, g, width = 8, height = 6, dpi = 300, units = "in", bg="white")

    df = rbind(df, data.frame(
      n.neighbors=n.neighbors,
      min.dist=min.dist,
      umap_png=umap_png,
      stringsAsFactors=FALSE
    ))
  }
}

write.csv(df, file = paste0("umap_png.csv"), row.names = FALSE)
