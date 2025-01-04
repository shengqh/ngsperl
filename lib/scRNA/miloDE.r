rm(list=ls()) 
sample_name='HiAdipose_vs_HiPBMC'
outFile='HiAdipose_vs_HiPBMC'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/seurat_sct_merge_dr0.5_3_choose/result/combined.final.rds'
parFile2=''
parFile3=''


setwd('/data/wanjalla_lab/projects/20230501_combined_scRNA_hg38/20241007_MiloDE/CMV_Bulk/result/HiAdipose_vs_HiPBMC')

### Parameter setting end ###

library(data.table)
library(ggplot2)
library(dichromat)
library(RColorBrewer)

options(future.globals.maxSize= 10779361280)

myoptions = fread(parSampleFile2, header=FALSE, data.table=FALSE)
myoptions = split(myoptions$V1, myoptions$V2)

ct = myoptions$celltype
cat("celltype =", ct, "\n")

reduction=myoptions$reduction
cat("reduction =", reduction, "\n")

cluster_id=myoptions$cluster_id
cat("cluster_id =", cluster_id, "\n")

optimized_neighbor_cells=as.numeric(myoptions$optimized_neighbor_cells)
cat("optimized_neighbor_cells =", optimized_neighbor_cells, "\n")

is_unix = .Platform$OS.type == "unix"

if(is_unix){
  library(BiocParallel)
  ncores = as.numeric(myoptions$ncores)
  mcparam = MulticoreParam(workers = ncores)
  register(mcparam)
}

library(Seurat)
library(miloDE)
library(SingleCellExperiment)

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

neightbour_file=paste0(prefix, ".neightbourhoods.rds")
if(!file.exists(neightbour_file)){
  sce_rds = paste0(prefix, ".sce.rds")
  if(file.exists(sce_rds)){
    comp_sce=readRDS(sce_rds)
  }else{
    cat("Reading object from file:", basename(parFile1), "\n")
    cur_obj=readRDS(parFile1)
    cur_obj=subset(cur_obj, sample %in% sample_groups$V1)
    cur_obj@meta.data$condition_id=unlist(sample_group_map[as.character(cur_obj$sample)])

    # Convert to SingleCellExperiment
    group_obj = subset(cur_obj, condition_id %in% groups$V1)
    rm(cur_obj)

    if(ct != "Bulk"){
      ct_obj=subset(group_obj, cell_type == ct)
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
  k_grid = seq(20, 100, 10)
  if(file.exists(size_rds)){
    stat_k = readRDS(size_rds)
  }else{
    cat("estimate_neighbourhood_sizes...\n")
    stat_k = estimate_neighbourhood_sizes(comp_sce, 
                                          k_grid = k_grid, 
                                          order = 2, 
                                          prop = 0.1 , 
                                          filtering = TRUE,
                                          cluster_id = cluster_id,
                                          reducedDim_name = reduction, 
                                          plot_stat = FALSE)
    saveRDS(stat_k, size_rds)
  }   

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

  #https://github.com/MarioniLab/miloDE/issues/29
  cat("assign_neighbourhoods\n")
  comp_sce = assign_neighbourhoods( comp_sce, 
                                    k = best_k,
                                    filtering = TRUE, 
                                    reducedDim_name = reduction)
  saveRDS(comp_sce, file=neightbour_file)

  unlink(sce_rds)
}else{
  comp_sce=readRDS(neightbour_file)
}

comp_sce$condition_id=factor(comp_sce$condition_id, levels=groups$V1)
comp_df = as.data.frame(table(comp_sce$condition_id)) |> 
  dplyr::mutate(Group=c("control", "sample")) |>
  dplyr::rename(Name=Var1, Cell=Freq) |>
  dplyr::select(Group, Name, Cell)
write.csv(comp_df, paste0(prefix, ".group_cell.csv"), row.names=FALSE)
print(comp_df)

de_file = paste0(prefix, ".miloDE.rds")
if(!file.exists(de_file)){
  cat("de_test_neighbourhoods\n")
  if(is_unix){
    de_stat = de_test_neighbourhoods( comp_sce, 
                                    sample_id = "sample", 
                                    design = ~condition_id, 
                                    covariates = c("condition_id"),
                                    BPPARAM = mcparam)
  }else{
    de_stat = de_test_neighbourhoods( comp_sce, 
                                    sample_id = "sample", 
                                    design = ~condition_id, 
                                    covariates = c("condition_id"))
  }
  saveRDS(de_stat, file=de_file)
}else{
  de_stat=readRDS(de_file)
}
