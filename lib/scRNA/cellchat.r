rm(list=ls()) 
sample_name='DMSO'
outFile='DMSO'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/20260410_silhouette_refine_clusters/Aorta_Progeria.final.obj.rds'
parFile2=''
parFile3=''


setwd('/nobackup/brown_lab/projects/20260324_Aorta_Progeria_scRNA_mm10/final_obj_cellchat/result/DMSO')

### Parameter setting end ###

source("scRNA_func.r")
#devtools::install_github("LTLA/BiocNeighbors")
#devtools::install_github("jinworks/CellChat")

library(data.table)
library(CellChat)
library(patchwork)
library(Seurat)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 100000000000)
options(parallelly.availableCores.logical=FALSE)

read_file_map<-function(file_list_path, sep="\t", header=FALSE) {
  tbl<-fread(file_list_path, header=header, data.table=FALSE, sep=sep)
  result<-split(tbl$V1, tbl$V2)
  return(result)
}

sample_list = read_file_map(parSampleFile1)
group_value = sample_list[[sample_name]]

myoptions = read_file_map(parSampleFile2)

future::plan("multisession", workers = as.numeric(myoptions$thread)) # do parallel

obj = readRDS(parFile1)
if(parFile2 != ""){
  obj@meta.data=readRDS(parFile2)
}

meta=obj@meta.data
group_meta=meta |> dplyr::filter(!!as.symbol(myoptions$group_column) == group_value)
if(myoptions$ignore_celltypes != ""){
  ignore_celltypes = trimws(unlist(strsplit(myoptions$ignore_celltypes, ",")))
  if(! all(ignore_celltypes %in% unique(group_meta[, myoptions$celltype_column]))) {
    missing_celltypes = setdiff(ignore_celltypes, unique(group_meta[, myoptions$celltype_column]))
    stop(paste0("The following cell types in ignore_celltypes are not found in the data: ", paste(missing_celltypes, collapse=", ")))  
  }
  cat(paste0("Ignoring cell types: ", paste(ignore_celltypes, collapse=", "), "\n"))
  group_meta = group_meta |> dplyr::filter(! (!!as.symbol(myoptions$celltype_column) %in% ignore_celltypes))
}
cells = rownames(group_meta)
obj = subset(obj, cells=cells)
obj = NormalizeData(obj)

data.input = GetAssayData(obj, assay = "RNA", layer = "data")

meta = obj@meta.data |>
  dplyr::rename(
    "labels"=myoptions$celltype_column,
    "samples"=myoptions$sample_column) |>
  dplyr::select(labels, samples) |>
  dplyr::mutate(
    labels = as.factor(as.character(labels)),
    samples = as.factor(as.character(samples))
 )

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

if(myoptions$CellChatDB == "mouse"){
  cellchat@DB <- CellChatDB.mouse
} else {
  cellchat@DB <- CellChatDB.human
}

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cat("identifyOverExpressedGenes ...\n")
cellchat <- identifyOverExpressedGenes(cellchat)

cat("identifyOverExpressedInteractions ...\n")
cellchat <- identifyOverExpressedInteractions(cellchat)

cat("computeCommunProb ...\n")
cellchat <- computeCommunProb(cellchat, type = "triMean")

cat("filterCommunication ...\n")
cellchat <- filterCommunication(cellchat, min.cells = 10)

cat("computeCommunProbPathway ...\n")
cellchat <- computeCommunProbPathway(cellchat)

cat("aggregateNet ...\n")
cellchat <- aggregateNet(cellchat)
    
cellchat_file=paste0(outFile, ".cellchat.rds")
saveRDS(cellchat, file = cellchat_file)

