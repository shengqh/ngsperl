rm(list=ls()) 
sample_name='CD_Met'
outFile='CD_Met'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1='/nobackup/h_cqs/ramirema/maureen_gannon/20240227_10940_snRNAseq_mmulatta_proteincoding_cellbender/P10940.RNA.ATAC.combined.finalnewcelltype.rds'
parFile2='/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/P10940.final.treatment.rds'
parFile3=''


setwd('/nobackup/h_cqs/maureen_gannon_projects/20240321_10940_snRNAseq_mmulatta_proteincoding_cellbender/20250428_cellchat/cellchat/result/CD_Met')

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
cells = rownames(obj@meta.data)[obj@meta.data[, myoptions$group_column] == group_value]
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

