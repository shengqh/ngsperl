rm(list=ls()) 
sample_name='ATC_vs_Normal'
outFile='ATC_vs_Normal'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parSampleFile4='fileList4.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20250403_fibroblasts_miloR_miloDE/milo_Bulk_2_miloR/result/ATC_vs_Normal')

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

source("reportFunctions.R")

myoptions = read_file_map(parSampleFile2, do_unlist=FALSE)

ct = get_option_value(myoptions, "celltype", required=TRUE)
cat("celltype =", ct, "\n")

#reduction for visualization only
visulization_reduction=get_option_value(myoptions, "visulization_reduction", required=TRUE)
cat("visulization_reduction =", visulization_reduction, "\n")

#neighbourhood_reduction for neighbourhood calculation, could be MNN for fastmnn (pca.corrected)
neighbourhood_reduction=get_option_value(myoptions, "neighbourhood_reduction", required=TRUE)
cat("neighbourhood_reduction =", neighbourhood_reduction, "\n")

#reductions would be uppercase in SingleCellExperiment
sce_visulization_reduction=toupper(visulization_reduction)
sce_neighbourhood_reduction=toupper(neighbourhood_reduction)

#annotation_column is used to annotate the neighbourhoods
annotation_column=get_option_value(myoptions, "annotation_column", required=FALSE)
cat("annotation_column =", annotation_column, "\n")

sub_annotation_column=get_option_value(myoptions, "sub_annotation_column", required=FALSE)
cat("sub_annotation_column =", sub_annotation_column, "\n")

SpatialFDR=as.numeric(get_option_value(myoptions, "SpatialFDR", required=FALSE, default=0.1))
cat("SpatialFDR =", SpatialFDR, "\n")

is_unix = .Platform$OS.type == "unix"

if(is_unix){
  library(BiocParallel)
  ncores = as.numeric(myoptions$ncores)
  mcparam = MulticoreParam(workers = ncores)
  register(mcparam)
}else{
  mcparam=NULL
}

pairs = fread(parSampleFile1, header=FALSE, data.table=FALSE)
cat("comparison =", sample_name, "\n")
groups = pairs |> dplyr::filter(V3==sample_name, V2=="groups")
cat("control groups =", groups$V1[1], "\n")
cat("case groups =", groups$V1[2], "\n")

prefix = paste0(sample_name, ".", celltype_to_filename(ct))
cat("Performing", prefix, "\n")

neighbourhoods_file=fread("fileList4.txt", header=FALSE, data.table=FALSE) |>
  dplyr::filter(V2==sample_name) |>
  dplyr::pull(V1)
cat("Reading object from file:", neighbourhoods_file, "\n")
comp_milo=readRDS(neighbourhoods_file)

comp_milo$milo_condition=factor(comp_milo$milo_condition, levels=groups$V1)

cat("countCells\n")
comp_milo <- countCells(comp_milo, 
                        meta.data = data.frame(colData(comp_milo)), 
                        samples="sample")

comp_design <- distinct(data.frame(colData(comp_milo))[,c("sample", "milo_condition")])
rownames(comp_design) <- comp_design$sample

cat("calcNhoodDistance\n")
comp_milo <- calcNhoodDistance( comp_milo, 
                                d = 30,
                                reduced.dim = sce_neighbourhood_reduction)

cat("testNhoods\n")
# Since we don't use contrast, we should not use design = ~ 0 + milo_condition which would cause all negative logFC
da_results <- testNhoods( comp_milo, 
                          design = ~ milo_condition, 
                          design.df = comp_design,
                          reduced.dim = sce_neighbourhood_reduction)

saveRDS(da_results, file=paste0(prefix, ".miloR_da.rds"))

cat("ploting\n")
#comp_milo <- buildNhoodGraph(comp_milo)

cat("\nplotNhoodGraphDA\n")
nh_graph_pl <- plotNhoodGraphDA(comp_milo, 
                                da_results, 
                                layout=sce_visulization_reduction, 
                                alpha=0.1) 
ggsave(paste0(prefix, ".nhood_graph.png"), nh_graph_pl, width = 8, height = 6, units = "in", dpi = 300, bg="white")


if(!is.null(annotation_column)){
  cat("\nplotNhoodGraphDA\n")
  da_results <- annotateNhoods( comp_milo, 
                                da_results, 
                                coldata_col = annotation_column)
  cluster_fraction=paste0(annotation_column, "_fraction")
  da_results[,paste0(annotation_column, "_final")] <- ifelse(da_results[, cluster_fraction] < 0.7, "Mixed", da_results[,annotation_column])
  g=ggplot(da_results, aes(!!sym(cluster_fraction))) + geom_histogram(bins=50)
  ggsave(paste0(prefix, ".annotation_fraction.png"), g, width = 4, height = 3, units = "in", dpi = 300, bg="white")
}

if(!is.null(sub_annotation_column)){
  da_results <- annotateNhoods( comp_milo, 
                                da_results, 
                                coldata_col = sub_annotation_column)
  sub_cluster_fraction=paste0(sub_annotation_column, "_fraction")
  da_results[,paste0(sub_annotation_column, "_final")] <- ifelse(da_results[,sub_cluster_fraction] < 0.7, "Mixed", da_results[,sub_annotation_column])
  g=ggplot(da_results, aes(!!sym(sub_cluster_fraction))) + geom_histogram(bins=50)
  ggsave(paste0(prefix, ".sub_annotation_fraction.png"), g, width = 4, height = 3, units = "in", dpi = 300, bg="white")
}

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
  if(!is.null(annotation_column)){
    g=plotDAbeeswarm(da_results, group.by =annotation_column, alpha=SpatialFDR) + theme(axis.title.y=element_blank())
    ggsave(paste0(prefix, ".annotation.DA_beeswarm.png"), g, width = 10, height = 5, units = "in", dpi = 300, bg="white")
  }

  if(!is.null(sub_annotation_column)){
    g=plotDAbeeswarm(da_results, group.by = sub_annotation_column, alpha=SpatialFDR) + theme(axis.title.y=element_blank())
    ggsave(paste0(prefix, ".sub_annotation.DA_beeswarm.png"), g, width = 10, height = 8, units = "in", dpi = 300, bg="white")
  }
}

