rm(list=ls()) 
sample_name='AS01'
outFile='AS01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20230512_Vandy_AS_scRNA_sct2_atac/archR/result/AS01')

### Parameter setting end ###

source("scRNA_func.r")

#https://www.archrproject.com/articles/Articles/tutorial.html
#https://github.com/GreenleafLab/ArchR_Website/blob/ui_updates/bookdown/21_MultiomicDataAnalysis.Rmd

library(ArchR)
set.seed(20230512)

myoptions<-read_file_map("fileList2.txt", do_unlist=FALSE)
minTSS=to_numeric(myoptions$minTSS, 4)
minFrags=to_numeric(myoptions$minFrags, 1000)
threads=to_numeric(myoptions$threads, 8)
genome=get_value(myoptions$genome, "hg38")
sampleCells=to_numeric(myoptions$sampleCells, 10000)
minDist=to_numeric(myoptions$minDist, 0.8)

lsi_resolution=to_numeric(myoptions$lsi_resolution, 0.2)
cluster_resolution=to_numeric(myoptions$cluster_resolution, 0.4)

addArchRThreads(threads = threads) 

inputFiles <- read_file_map("fileList1.txt", do_unlist=TRUE)
inputFiles

addArchRGenome(genome)

#Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = minTSS, #Dont set this too high because you can always increase later
  minFrags = minFrags, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = ".",
  copyArrows = FALSE #
  #copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#Filtering doublets
proj <- addDoubletScores(proj, useMatrix = "TileMatrix", force = TRUE)
proj <- filterDoublets(proj)

#getAvailableMatrices(proj)

#Dimensionality Reduction and Clustering

#ArchR implements an iterative LSI dimensionality reduction
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = lsi_resolution, 
    sampleCells = sampleCells,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

proj <- addClusters(input = proj, reducedDims = "LSI_ATAC")

#Visualizing in a 2D UMAP Embedding
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = minDist, force = TRUE)

proj <- addClusters(proj, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = cluster_resolution, force = TRUE)

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters_ATAC", embedding = "UMAP_ATAC") + theme(aspect.ratio=1)
png(paste0(outFile, ".cluster.umap.png"), width=3000, height=2500, res=300)
print(p1)
dev.off()

if(length(inputFiles) > 0){
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP_ATAC") + theme(aspect.ratio=1)
  png(paste0(outFile, ".sample.umap.all.png"), width=3000, height=2500, res=300)
  print(p1)
  dev.off()
}

#saveArchRProject(proj, outFile)

#Assigning Clusters with Gene Scores





