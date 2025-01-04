rm(list=ls()) 
sample_name='HiCMV_vs_LoCMV'
outFile='HiCMV_vs_LoCMV'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/seurat_fastmnn_dr0.5_3_choose/result/combined.final.rds'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38_fastmnn/20241230_miloR/CMV_Bulk/result/HiCMV_vs_LoCMV')

### Parameter setting end ###

library(data.table)
library(ggplot2)
library(dichromat)
library(RColorBrewer)
library(Seurat)
library(SingleCellExperiment)
library(miloR)
library(miloDE)

source('scRNA_func.r')

options(future.globals.maxSize= 10779361280)

myoptions = fread(parSampleFile2, header=FALSE, data.table=FALSE)
myoptions = split(myoptions$V1, myoptions$V2)

ct = myoptions$celltype
cat("celltype =", ct, "\n")

reduction=myoptions$reduction
cat("reduction =", reduction, "\n")

pca=myoptions$pca
cat("pca =", pca, "\n")

cluster_id=myoptions$cluster_id
cat("cluster_id =", cluster_id, "\n")

sub_cluster_id=myoptions$sub_cluster_id
cat("sub_cluster_id =", sub_cluster_id, "\n")

SpatialFDR=as.numeric(myoptions$SpatialFDR)
cat("SpatialFDR =", SpatialFDR, "\n")

optimized_neighbor_cells=as.numeric(myoptions$optimized_neighbor_cells)
cat("optimized_neighbor_cells =", optimized_neighbor_cells, "\n")

merge_microphages_DC=ifelse("merge_macrophages_DC" %in% names(myoptions), myoptions$merge_macrophages_DC=="1", FALSE)
cat("merge_microphages_DC =", merge_microphages_DC, "\n")

is_unix = .Platform$OS.type == "unix"

if(is_unix){
  library(BiocParallel)
  ncores = as.numeric(myoptions$ncores)
  mcparam = MulticoreParam(workers = ncores)
  register(mcparam)
}

pairs = fread(parSampleFile1, header=FALSE, data.table=FALSE)
cat("comparison =", sample_name, "\n")
groups = pairs |> dplyr::filter(V3==sample_name, V2=="groups")
cat("control groups =", groups$V1[1], "\n")
cat("sample groups =", groups$V1[2], "\n")

sample_groups = fread(parSampleFile3, header=FALSE, data.table=FALSE) 
sample_groups = sample_groups |> dplyr::filter(V2 %in% groups$V1)
sample_group_map=split(sample_groups$V2, sample_groups$V1)

prefix = gsub(" ", "_", paste0(sample_name, ".", ct))
cat("Performing", prefix, "\n")

sce_obj_file=paste0(prefix, ".sce.rds")
if(!file.exists(sce_obj_file)){
  cat("Reading object from file:", basename(parFile1), "\n")
  cur_obj=readRDS(parFile1)
  cur_obj=subset(cur_obj, sample %in% sample_groups$V1)
  cur_obj@meta.data$condition_id=unlist(sample_group_map[as.character(cur_obj$sample)])

  group_obj = subset(cur_obj, condition_id %in% groups$V1)
  group_obj@meta.data$condition_id=factor(group_obj@meta.data$condition_id, levels=groups$V1)
  rm(cur_obj)

  if(ct != "Bulk"){
    ct_obj=subset(group_obj, !!sym(cluster_id) == ct)
  }else{
    ct_obj=group_obj
  }
  rm(group_obj)

  ct_obj@meta.data$SummaryLayer=as.character(ct_obj@meta.data[,cluster_id])
  if(merge_microphages_DC){
    ct_obj@meta.data$SummaryLayer[ct_obj@meta.data$SummaryLayer=="Macrophages"]="Macrophages/DC"
    ct_obj@meta.data$SummaryLayer[ct_obj@meta.data$SummaryLayer=="Dendritic cells"]="Macrophages/DC"
  }

  discard_ct=table(ct_obj$SummaryLayer) |> 
    as.data.frame() |> 
    dplyr::filter(Freq < 10) |> 
    dplyr::pull(Var1)

  if(length(discard_ct) > 0) {
    cat("Discarding", paste0(discard_ct, collapse=","), "with less than 10 cells\n")
    ct_obj=subset(ct_obj, SummaryLayer %in% discard_ct, invert=TRUE)
  }
  ct_obj@meta.data$SummaryLayer=factor_by_count(ct_obj@meta.data$SummaryLayer)

  g<-get_dim_plot_labelby(ct_obj, label.by="SummaryLayer", reduction=reduction, ggplot_default_colors=TRUE) + 
    theme(legend.title=element_blank(),
          plot.title=element_blank(),
          axis.ticks=element_blank()) 
  ggsave(paste0(prefix, ".UMAP.png"), g, width = 6.4, height = 4.4, units = "in", dpi = 300, bg="white")

  g<-DimPlot(ct_obj, group.by="condition_id", label=FALSE, reduction=reduction) + 
    facet_wrap(~condition_id, ncol=1) + 
    theme(legend.position="none",
          plot.title=element_blank(),
          axis.line=element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(),
          strip.background=element_blank(),
          aspect.ratio=1)
  ggsave(paste0(prefix, ".UMAP.condition.png"), g, width = 2, height = 4.2, units = "in", dpi = 300, bg="white")

  # g=DimPlot_scCustom(seurat_object = ct_obj, group.by = "condition_id")
  # ggsave(paste0(prefix, ".UMAP.condition_id.png"), g, width = 6, height = 4, units = "in", dpi = 300, bg="white")

  cat("as.SingleCellExperiment\n")
  comp_sce=as.SingleCellExperiment(ct_obj)
  rm(ct_obj)
  cat("Saving object to file:", sce_obj_file, "\n")
  saveRDS(comp_sce, file=sce_obj_file)
}else{
  cat("Reading object from file:", sce_obj_file, "\n")
  comp_sce=readRDS(sce_obj_file)
}

#reductions would be uppercase in SingleCellExperiment
sce_reduction=toupper(reduction)
sce_pca=toupper(pca)

comp_df = as.data.frame(table(comp_sce$condition_id)) |> 
  dplyr::mutate(Group=c("control", "sample")) |>
  dplyr::rename(Name=Var1, Cell=Freq) |>
  dplyr::select(Group, Name, Cell)
write.csv(comp_df, paste0(prefix, ".group_cell.csv"), row.names=FALSE)
print(comp_df)

k_grid = seq(20, 100, 10)
cat("estimate_neighbourhood_sizes...\n")
stat_k = estimate_neighbourhood_sizes(comp_sce, 
                                      k_grid = k_grid, 
                                      order = 2, 
                                      prop = 0.1 , 
                                      filtering = TRUE,
                                      cluster_id = cluster_id,
                                      reducedDim_name = sce_pca, 
                                      plot_stat = FALSE)
saveRDS(stat_k, paste0(prefix, ".estimate_sizes.rds"))

best_k = stat_k |> 
  dplyr::mutate(distance=abs(med-optimized_neighbor_cells)) |>
  dplyr::slice(which.min(distance)) |>
  dplyr::pull(k) |> 
  as.character() |> 
  as.numeric()

cat("Pick the best k as", best_k, "for median of number of cells closest to", optimized_neighbor_cells, "\n")  
writeLines(as.character(best_k), paste0(prefix, ".best_k.txt"))

size_png = paste0(prefix, '.estimate_sizes.png')
g = ggplot(stat_k, aes(k)) + 
    geom_boxplot(aes(ymin = min, lower = q25, middle = med, upper = q75, ymax = max, fill = k), stat = "identity") + 
    theme_bw() + 
    scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(length(k_grid))) + 
    labs(y = "Neighbourhood size")
ggsave(size_png, g, width = 6, height = 4, units = "in", dpi = 300, bg="white")

cat("assign_neighbourhoods\n")
comp_milo = assign_neighbourhoods( comp_sce, 
                                  k = best_k,
                                  filtering = TRUE, 
                                  reducedDim_name = sce_pca)

nhood_png=paste0(prefix, ".nhood.png")
g=plotNhoodSizeHist(comp_milo)
ggsave(nhood_png, g, width = 4, height = 3, units = "in", dpi = 300, bg="white")

comp_milo <- countCells(comp_milo, meta.data = as.data.frame(colData(comp_milo)), sample="sample")

nhood_count_csv=paste0(prefix, ".nhood_count.csv")
write.csv(nhoodCounts(comp_milo), nhood_count_csv)

comp_design <- distinct(data.frame(colData(comp_milo))[,c("sample", "condition_id")])
rownames(comp_design) <- comp_design$sample

cat("calcNhoodDistance\n")
comp_milo <- calcNhoodDistance(comp_milo, d=30, reduced.dim = sce_pca)

saveRDS(comp_milo, file=paste0(prefix, ".miloR.rds"))

cat("testNhoods\n")
da_results <- testNhoods(comp_milo, design = ~ condition_id, design.df = comp_design)

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

cat("ploting\n")
comp_milo <- buildNhoodGraph(comp_milo)

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(comp_milo, da_results, layout=sce_reduction, alpha=0.1) 
ggsave(paste0(prefix, ".nhood_graph.png"), nh_graph_pl, width = 9, height = 7, units = "in", dpi = 300, bg="white")

da_results <- annotateNhoods(comp_milo, da_results, coldata_col = "SummaryLayer")
da_results[,"SummaryLayer_final"] <- ifelse(da_results[,"SummaryLayer_fraction"] < 0.7, "Mixed", da_results[,"SummaryLayer"])

g=ggplot(da_results, aes(SummaryLayer_fraction)) + geom_histogram(bins=50)
ggsave(paste0(prefix, ".cluster_fraction.png"), g, width = 4, height = 3, units = "in", dpi = 300, bg="white")

sub_cluster_fraction=paste0(sub_cluster_id, "_fraction")
da_results <- annotateNhoods(comp_milo, da_results, coldata_col = sub_cluster_id)

da_results[,paste0(sub_cluster_id, "_final")] <- ifelse(da_results[,sub_cluster_fraction] < 0.7, "Mixed", da_results[,sub_cluster_id])
g=ggplot(da_results, aes(!!sym(sub_cluster_fraction))) + geom_histogram(bins=50)
ggsave(paste0(prefix, ".sub_cluster_fraction.png"), g, width = 4, height = 3, units = "in", dpi = 300, bg="white")

da_results = da_results |> dplyr::arrange(SpatialFDR)
write.csv(da_results, paste0(prefix, ".da_results.csv"), row.names=FALSE)

g=ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) + theme_bw()
ggsave(paste0(prefix, ".PValue.png"), g, width = 4, height = 3, units = "in", dpi = 300, bg="white")

g=ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = -log10(SpatialFDR)) +
  theme_bw()
ggsave(paste0(prefix, ".volcano.png"), g, width = 6, height = 5, units = "in", dpi = 300, bg="white")

if(any(da_results$SpatialFDR <= SpatialFDR)){
  g=plotDAbeeswarm(da_results, group.by ="SummaryLayer", alpha=SpatialFDR) + theme(axis.title.y=element_blank())
  ggsave(paste0(prefix, ".cluster.DA_beeswarm.png"), g, width = 10, height = 5, units = "in", dpi = 300, bg="white")

  g=plotDAbeeswarm(da_results, group.by = sub_cluster_id, alpha=SpatialFDR) + theme(axis.title.y=element_blank())
  ggsave(paste0(prefix, ".sub_cluster.DA_beeswarm.png"), g, width = 10, height = 8, units = "in", dpi = 300, bg="white")
}

unlink(sce_obj_file)
