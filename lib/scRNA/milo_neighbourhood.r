rm(list=ls()) 
sample_name='ATC_vs_Normal'
outFile='ATC_vs_Normal'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/data/h_vivian_weiss/Thyroid_scRNA_seq_atlas_files/24-0819_Atlas/24-1114_Atlas_Files/24-0821_Fibroblast_Subclustering_FastMNN_3000.RDS'
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/shengq2/20250403_fibroblasts_miloR_miloDE/milo_Bulk_1_neighbourhood/result/ATC_vs_Normal')

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

#https://rawcdn.githack.com/MarioniLab/miloDE_tutorials/3d3781237011695f802dc1c0f0193bea12a108de/miloDE__mouse_embryo.html
#get estimated neighbourhood sizes which are close to optimized_neighbour_cells (400 cells in toturial)
optimized_neighbour_cells=as.numeric(get_option_value(myoptions, "optimized_neighbour_cells", required=FALSE, default=400))
cat("optimized_neighbour_cells =", optimized_neighbour_cells, "\n")

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

prefix = gsub(" ", "_", paste0(sample_name, ".", ct))
cat("Performing", prefix, "\n")

neighbour_file=paste0(prefix, ".milo.neighbourhoods.rds")
if (file.exists(neighbour_file)){
  cat("Reading object from file:", basename(neighbour_file), "\n")
  comp_milo=readRDS(neighbour_file)
}else{
  sce_rds = paste0(prefix, ".sce.rds")
  if(file.exists(sce_rds)){
    comp_sce=readRDS(sce_rds)
  }else{
    cat("Reading object from file:", basename(parFile1), "\n")
    cur_obj=readRDS(parFile1)

    if(!"sample" %in% colnames(cur_obj@meta.data)){
      cur_obj@meta.data$sample=cur_obj@meta.data$orig.ident
    }

    if(is.null(myoptions$condition_column)){
      sample_groups = fread(parSampleFile3, header=FALSE, data.table=FALSE) 
      sample_groups = sample_groups |> dplyr::filter(V2 %in% groups$V1)
    }else{
      sample_groups = cur_obj@meta.data[,c("orig.ident", myoptions$condition_column)] |> unique() |>
        dplyr::rename(V1=1, V2=2) |>
        dplyr::mutate(V1=as.character(V1),
                      V2=as.character(V2))
      sample_groups = sample_groups |> dplyr::filter(V2 %in% groups$V1)
    }
    sample_group_map=split(sample_groups$V2, sample_groups$V1)

    cur_obj=subset(cur_obj, sample %in% sample_groups$V1)
    cur_obj@meta.data$milo_condition=unlist(sample_group_map[as.character(cur_obj$sample)])

    # Convert to SingleCellExperiment
    group_obj = subset(cur_obj, milo_condition %in% groups$V1)
    rm(cur_obj)

    if(ct != "Bulk"){
      ct_obj=subset(group_obj, !!sym(annotation_column) == ct)
    }else{
      ct_obj=group_obj
    }

    cat("as.SingleCellExperiment\n")
    comp_sce=as.SingleCellExperiment(ct_obj)
    saveRDS(comp_sce, file=sce_rds)

    rm(group_obj)
    rm(ct_obj)
  }

  size_rds =  paste0(prefix, ".estimate_sizes.rds")
  k_grid = seq(10, 50, 10)
  if(file.exists(size_rds)){
    stat_k = readRDS(size_rds)
  }else{
    #since cluster_id is an optinal parameter in estimate_neighbourhood_sizes but not in assign_neighbourhoods, 
    #I am not sure if we should use it. So I removed it from estimate_neighbourhood_sizes
    cat("estimate_neighbourhood_sizes...\n")
    stat_k = estimate_neighbourhood_sizes(comp_sce, 
                                          reducedDim_name = sce_neighbourhood_reduction, 
                                          k_grid = k_grid, 
                                          order = 2, 
                                          prop = 0.1 , 
                                          filtering = TRUE,
                                          plot_stat = FALSE)
    saveRDS(stat_k, size_rds)
  }   

  best_k = stat_k |> 
    dplyr::mutate(distance=abs(med-optimized_neighbour_cells)) |>
    dplyr::slice(which.min(distance)) |>
    dplyr::pull(k) |> 
    as.character() |> 
    as.numeric()

  cat("Pick the best k as", best_k, "for median of number of cells closest to", optimized_neighbour_cells, "\n")  
  writeLines(as.character(best_k), paste0(prefix, ".best_k.txt"))

  g = ggplot(stat_k, aes(k)) + 
      geom_boxplot(aes(ymin = min, lower = q25, middle = med, upper = q75, ymax = max, fill = k), stat = "identity") + 
      theme_bw() + 
      scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(length(k_grid))) + 
      labs(y = "Neighbourhood size")
  ggsave(paste0(prefix, '.estimate_sizes.png'), g, width = 6, height = 4, units = "in", dpi = 300, bg="white")

  #https://github.com/MarioniLab/miloDE/issues/29
  cat("assign_neighbourhoods\n")
  #It calls miloR::buildGraph internally
  comp_milo = assign_neighbourhoods(comp_sce, 
                                    reducedDim_name = sce_neighbourhood_reduction,
                                    k = best_k,
                                    filtering = TRUE )

  cat("save neighbourhoods object\n")
  saveRDS(comp_milo, file=neighbour_file)
  unlink(sce_rds)
}

cdf = data.frame(Group=c("control", "case"), Name=groups$V1)
comp_df = as.data.frame(table(comp_milo$milo_condition)) |> 
  dplyr::rename(Name=Var1, Cell=Freq) |>
  dplyr::left_join(cdf, by="Name") |>
  dplyr::select(Group, Name, Cell)
write.csv(comp_df, paste0(prefix, ".group_cell.csv"), row.names=FALSE)

cat("plotNhoodSizeHist\n")
g=plotNhoodSizeHist(comp_milo)
ggsave(paste0(prefix, ".nhood_hist.png"), g, width = 5, height = 4, units = "in", dpi = 300, bg="white")

if(!is.null(annotation_column)) {
  cat("annotateNhoods\n")
  nhood_stat_ct = get_nhood_stat_ct(comp_milo, annotation_column)
  write.csv(nhood_stat_ct, paste0(prefix, ".nhoods_annotation.csv"), row.names=FALSE)

  cat("plot_milo_by_single_metric\n")
  p = plot_milo_by_single_metric( comp_milo, 
                                  nhood_stat_ct, 
                                  colour_by = annotation_column , 
                                  layout = sce_visulization_reduction , 
                                  size_range = c(1.5,3) , 
                                  edge_width = c(0.2,0.5)) +
      theme(aspect.ratio=1)

  nhoods_width=12
  nhoods_height=8

  ggsave(paste0(prefix, ".nhoods_annotation.png"), p, width = nhoods_width, height = nhoods_height, units = "in", dpi = 300, bg="white")
}

cat("done\n")
