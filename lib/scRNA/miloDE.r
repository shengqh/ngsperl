rm(list=ls()) 
sample_name='ATC_vs_Normal'
outFile='ATC_vs_Normal'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/nobackup/h_vivian_weiss_lab/shengq2/20250401_atlas_miloR_miloDE/Bulk_1_milo_neighbourhood/result/ATC_vs_Normal/ATC_vs_Normal.Bulk.milo.neighbourhoods.rds'
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20250401_atlas_miloR_miloDE/Bulk_3_miloDE/result/ATC_vs_Normal')

### Parameter setting end ###

library(data.table)
library(ggplot2)
library(dichromat)
library(RColorBrewer)
library(Seurat)
library(miloDE)
library(miloR)
library(SingleCellExperiment)
library(viridis)

source("scRNA_func.r")

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

is_unix = .Platform$OS.type == "unix"
if(is_unix){
  library(BiocParallel)
  ncores = as.numeric(myoptions$ncores)
  mcparam = MulticoreParam(workers = ncores)
  register(mcparam)
}else{
  mcparam=NULL
}

set.seed(20250402)

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

if(!is.null(myoptions$filter_by_AUC) & myoptions$filter_by_AUC=="1"){
  stat_acu_rds=paste0(prefix, ".stat_auc.rds")
  if(file.exists(stat_acu_rds)){
    stat_auc = readRDS(stat_acu_rds)
  }else{
    cat("calc_AUC_per_neighbourhood\n")
    stat_auc = suppressWarnings(calc_AUC_per_neighbourhood( comp_milo, 
                                                            sample_id = "sample" , 
                                                            condition_id = "milo_condition", 
                                                            min_n_cells_per_sample = 1, 
                                                            BPPARAM = mcparam))
    saveRDS(stat_auc, file=stat_acu_rds)
  }

  g = plot_milo_by_single_metric( comp_milo, 
                                  stat_auc, 
                                  colour_by = "auc" , 
                                  layout = sce_visulization_reduction , 
                                  size_range = c(1.5,3) , 
                                  edge_width = c(0.2,0.5)) + 
                                  scale_fill_viridis(name = "AUC")

  ggsave(paste0(prefix, '.auc.png'), g, width = 6, height = 4, units = "in", dpi = 300, bg="white")
  # https://rawcdn.githack.com/MarioniLab/miloDE_tutorials/3d3781237011695f802dc1c0f0193bea12a108de/miloDE__mouse_embryo.html
  subset_nhoods=stat_auc$Nhood[!is.na(stat_auc$auc) & stat_auc$auc >= 0.5]
}else{
  subset_nhoods=NULL
}

de_file = paste0(prefix, ".miloDE.rds")
if(!file.exists(de_file)){
  cat("de_test_neighbourhoods\n")
  de_stat = de_test_neighbourhoods( comp_milo, 
                                    sample_id = "sample", 
                                    design = ~milo_condition, 
                                    covariates = c("milo_condition"),
                                    BPPARAM = mcparam)
  saveRDS(de_stat, file=de_file)
}else{
  de_stat=readRDS(de_file)
}

valid_de=de_stat |>
  dplyr::filter(!is.na(pval))

logFC_criteria=log2(1.5)

sig_de = valid_de |>
  dplyr::filter(pval_corrected_across_genes < 0.05) |>
  dplyr::filter(abs(logFC) > logFC_criteria) |>
  dplyr::select(-test_performed) |>
  dplyr::arrange(Nhood, pval)

if(!is.null(annotation_column)) {
  cat("annotateNhoods\n")
  nhood_stat_ct = get_nhood_stat_ct(comp_milo, annotation_column)

  sig_de_df = merge(sig_de, 
                    nhood_stat_ct |> dplyr::select(-Nhood_center), 
                    by="Nhood") 
}else{
  sig_de_df = sig_de
}

write.csv(sig_de_df, paste0(prefix, ".miloDE.significant.csv"), row.names=FALSE)

sig_genes = sig_de |> 
  dplyr::distinct(gene) |>
  dplyr::pull(gene) |>
  sort()

writeLines(sig_genes, paste0(prefix, ".miloDE.significant_genes.txt"))

cat("draw volcano plot...\n")
volcano_width=7
volcano_height=7

library(EnhancedVolcano)

p<-EnhancedVolcano(valid_de,
  lab = valid_de$gene,
  x = 'logFC',
  y = 'pval',
  title = sample_name,
  pCutoff = 0.05,
  pCutoffCol = 'pval_corrected_across_genes',
  FCcutoff = logFC_criteria,
  pointSize = 3.0,
  labSize = 6.0,
  colAlpha = 1,
  subtitle = NULL)

ggsave(paste0(prefix, ".miloDE.volcano_all.png"), p, width = volcano_width, height = volcano_height, units = "in", dpi = 300, bg="white")

unique_valid_de = valid_de |>
  dplyr::arrange(gene, pval) |>
  dplyr::distinct(gene, .keep_all=TRUE)

p<-EnhancedVolcano(unique_valid_de,
  lab = unique_valid_de$gene,
  x = 'logFC',
  y = 'pval',
  title = sample_name,
  pCutoff = 0.05,
  pCutoffCol = 'pval_corrected_across_genes',
  FCcutoff = logFC_criteria,
  pointSize = 3.0,
  labSize = 6.0,
  colAlpha = 1,
  subtitle = NULL)

ggsave(paste0(prefix, ".miloDE.volcano_unique.png"), p, width = volcano_width, height = volcano_height, units = "in", dpi = 300, bg="white")

stat_de_magnitude = rank_neighbourhoods_by_DE_magnitude(de_stat)

p1 = plot_milo_by_single_metric(comp_milo, 
                                stat_de_magnitude, 
                                colour_by = "n_DE_genes" , 
                                layout = sce_visulization_reduction , 
                                size_range = c(1.5,3) , 
                                edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# DE genes") +
  theme(aspect.ratio=1)
#> Warning in plot_milo_by_single_metric(comp_milo, stat_de_magnitude, colour_by =
#> "n_DE_genes", : Coercing layout to matrix format
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.
p2 = plot_milo_by_single_metric(comp_milo, 
                                stat_de_magnitude, 
                                colour_by = "n_specific_DE_genes" , 
                                layout = sce_visulization_reduction , 
                                size_range = c(1.5,3) , 
                                edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# specific\nDE genes" , option = "inferno")+
  theme(aspect.ratio=1)
#> Warning in plot_milo_by_single_metric(comp_milo, stat_de_magnitude, colour_by =
#> "n_specific_DE_genes", : Coercing layout to matrix format
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.
p = ggarrange(p1,p2)

num_gene_png=paste0(prefix, ".miloDE.num_gene.png")
ggsave(num_gene_png, p, width = 12, height = 5, units = "in", dpi = 300, bg="white")

cat("done\n")
