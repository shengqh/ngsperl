rm(list=ls()) 
sample_name='S03_ClassPTC_BRAF'
outFile='S03_ClassPTC_BRAF'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_vivian_weiss_lab/12904_RB_VisiumHD/20260212_12904_VisiumHD_cellsegment/MEcell_cluster_Leiden_report_subcluster/result/S03_ClassPTC_BRAF')

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

cluster_map<-fread("fileList3.txt", sep="\t", header=FALSE, fill=TRUE)
sub_resolutions=cluster_map |> dplyr::filter(V2=="resolutions") |> dplyr::pull(V1) |> as.numeric()

cluster_map=cluster_map |> dplyr::filter(V2!="resolutions")

cluster_map=cluster_map |> dplyr::filter(V3==sample_name)
resolution=cluster_map |> dplyr::filter(V2=="resolution") |> dplyr::pull(V1) |> as.numeric()
rename_map=cluster_map |> dplyr::filter(V2!="resolution")

log_info(paste0("Load Seurat object from ", rds_file, " ..."))
obj = readRDS(rds_file)
obj = UpdateSeuratObject(obj)

resolution_col=paste0("MEcell_res.", resolution)

stopifnot(assay %in% names(obj@assays))

DefaultAssay(obj) <- assay
log_info(paste0("There are ", ncol(obj), " cells/spots in the Seurat object."))

if(any(duplicated(rename_map$V1))){
  dup_clusters=rename_map$V1[duplicated(rename_map$V1)]
  stop(paste0("There are duplicated cluster numbers in rename_map: ", paste0(dup_clusters, collapse=", ")))
}

missed=setdiff(rename_map$V1, c(0, unique(as.numeric(obj@meta.data[,resolution_col]))))
if(length(missed) > 0){
  stop(paste0("Not all clusters in rename_map are present in the Seurat object metadata column ", resolution_col, ": ", paste0(missed, collapse=", ")))
}

algorithm=as.numeric(myoptions$cluster_algorithm)
algorithm_name=myoptions$cluster_algorithm_name

log_info(paste0("Perform MEcell sub-clustering with resolutions ", paste0(sub_resolutions, collapse=", "), " using algorithm ", algorithm_name))

celltypes=sort(unique(rename_map$V2))
celltype=celltypes[2]

res_tbl=NULL

obj@meta.data$previous_clusters=obj@meta.data[,resolution_col]
obj@meta.data$previous_celltype="Unknown"
obj@meta.data = obj@meta.data |> dplyr::select(!starts_with("MEcell_"))
for(celltype in celltypes){
  ct_name = gsub("[\\/\\s]+", "", celltype)
  ct_clusters=rename_map |> dplyr::filter(V2==celltype) |> dplyr::pull(V1)
  if(length(ct_clusters) == 1 && ct_clusters[1] == 0) {
    all_others=rename_map |> dplyr::filter(V2!=celltype) |> dplyr::pull(V1)
    ct_clusters=sort(setdiff(unique(as.numeric(obj@meta.data$previous_clusters)), all_others))
  }
  ct_cells=rownames(obj@meta.data[obj@meta.data$previous_clusters %in% ct_clusters,,drop=FALSE])

  obj@meta.data[ct_cells, "previous_celltype"] = celltype

  log_info(paste0(celltype, ": " , paste0(ct_clusters, collapse=", "), ": ", length(ct_cells), " cells, perform subclustering ..."))

  sub_obj=subset(obj, cells=ct_cells)

  sub_obj <- FindClusters(sub_obj, resolution = sub_resolutions, graph.name="MEcell", algorithm=algorithm, random.seed = 20251026)

  res_rds = paste0(sample_name, ".", ct_name, ".", algorithm_name, ".meta.rds")

  log_info(paste0("Save Seurat object metadata with MEcell clustering results to ", res_rds))
  saveRDS(sub_obj@meta.data, res_rds)

  cur_res_tbl=data.frame(celltype=celltype, cluster=paste0(ct_clusters, collapse=","), sub_resolutions=paste0(sub_resolutions, collapse=","), n_cells=length(ct_cells), meta_file=res_rds, stringsAsFactors=FALSE)
  if(is.null(res_tbl)){
    res_tbl=cur_res_tbl
  } else {
    res_tbl=rbind(res_tbl, cur_res_tbl)
  }
}

saveRDS(obj@meta.data, paste0(sample_name, ".meta.rds"))

write.csv(res_tbl, paste0(sample_name, ".", algorithm_name, ".subcluster.files.csv"), row.names=FALSE)

log_info("Done.")

