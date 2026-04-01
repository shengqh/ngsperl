rm(list=ls()) 
sample_name='S19_DSVPTC_CCND1RET'
outFile='S19_DSVPTC_CCND1RET'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/12904_RB_VisiumHD/20260212_12904_VisiumHD_cellsegment/MEcell_cluster_Leiden_report_subcluster_choose/result/S19_DSVPTC_CCND1RET')

### Parameter setting end ###

source("reportFunctions.R")
source("scRNA_func.r")
library('sf')
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(logger)

log_appender(appender_tee(file(paste0(sample_name, ".log"), open = "w")))

options(future.globals.maxSize = 600000 * 1024^2)
options(future.seed = TRUE)

myoptions<-read_file_map("fileList2.txt", do_unlist=FALSE)
email=myoptions$email
affiliation=get_affiliation(myoptions)
assay=myoptions$assay
species=myoptions$species

file_map<-read_file_map("fileList1.txt", do_unlist=FALSE)
sample_name=names(file_map)[1]
rds_file=file_map[[sample_name]]

meta_rds=paste0(sample_name, ".all.meta.rds")
# If we change the parameters, we will need to re-run the code. 
# So we will not skip the code even if the meta_rds file exists. 
# if(file.exists(meta_rds)){
#   quit(save="no")
# }

cluster_map<-fread("fileList3.txt", sep="\t", header=FALSE, fill=TRUE, data.table=FALSE) |>
  dplyr::filter(V4 == sample_name) |>
  dplyr::select(-V4) |>
  dplyr::mutate(V1=as.character(V1))

meta_file_list=fread('fileList4.txt', sep="\t", header=FALSE, fill=TRUE, data.table=FALSE) |> 
  dplyr::filter(V2==sample_name) |> dplyr::pull(V1)
meta_files=fread(meta_file_list, header=TRUE, data.table=FALSE)
meta_files$meta_file=paste0(dirname(meta_file_list), "/", meta_files$meta_file)

meta_cts=sort(unique(meta_files$celltype))
cat(paste0("Found cell types in subclustering: ", paste(meta_cts, collapse=", "), "\n"))

choose_cts=sort(unique(cluster_map$V3))
cat(paste0("Found cell types in cluster choose definition: ", paste(choose_cts, collapse=", "), "\n"))

if(!all(choose_cts %in% meta_cts)){
  stop(paste0("Not all cell types in cluster choose are found in the subclustering files. Missing: ", paste(choose_cts[!choose_cts %in% meta_cts], collapse=", ")))
}

if(!all(meta_cts %in% choose_cts)){
  stop(paste0("Not all cell types in subclustering files are found in the cluster choose. Missing: ", paste(meta_cts[!meta_cts %in% choose_cts], collapse=", ")))
}

log_info(paste0("Load Seurat object from ", rds_file, " ..."))
obj = readRDS(rds_file)
obj = UpdateSeuratObject(obj)

meta=obj@meta.data
meta$cell_type="DELETE"

ct=choose_cts[1]
for(ct in choose_cts) {
  cat(paste0("Processing cell type ", ct, " ...\n"))
  cur_clusters=cluster_map |> dplyr::filter(V3==ct) 
  cur_res=cur_clusters |> dplyr::filter(V2=="resolution") |> dplyr::pull(V1) |> as.numeric()
  cur_clusters=cur_clusters |> dplyr::filter(V2!="resolution")
  ct_meta_file=meta_files |> dplyr::filter(celltype==ct) |> dplyr::pull(meta_file)
  cat(paste0("Reading meta file ", ct_meta_file, " ...\n"))
  ct_meta=readRDS(ct_meta_file)
  ct_meta$cell_type="DELETE"

  cur_col=paste0("MEcell_res.", cur_res)
  if(!cur_col %in% colnames(ct_meta)){
    stop(paste0("Column ", cur_col, " not found in metadata file ", ct_meta_file))
  }

  if("-1" %in% cur_clusters$V1){
    default_ct=cur_clusters |> dplyr::filter(V1 == -1) |> dplyr::pull(V2) |> unique()
    ct_meta$cell_type=default_ct
    cur_clusters=cur_clusters |> dplyr::filter(V1 != -1)
  }

  if(nrow(cur_clusters) > 0) {
    i=1
    for(i in 1:nrow(cur_clusters)){
      cluster_num=cur_clusters$V1[i]
      cluster_num=strsplit(cluster_num, ",")[[1]]
      cluster_ct=cur_clusters$V2[i]
      ct_meta$cell_type[ct_meta[,cur_col] %in% cluster_num]=cluster_ct
    }
  }

  meta[rownames(ct_meta), "cell_type"]=ct_meta$cell_type
}

meta$cell_type=factor_by_count(meta$cell_type)

saveRDS(meta, meta_rds)

log_info("Done.")

