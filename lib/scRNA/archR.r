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

library(ArchR)
set.seed(20230512)

myoptions<-read_file_map("fileList2.txt", do_unlist=FALSE)
minTSS=to_numeric(myoptions$minTSS, 4)
minFrags=to_numeric(myoptions$minFrags, 1000)
threads=to_numeric(myoptions$threads, 8)
genome=get_value(myoptions$genome, "hg38")

addArchRThreads(threads = threads) 

inputFiles <- read_file_map("fileList1.txt", do_unlist=TRUE)
inputFiles

addArchRGenome(genome)

#Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#Inferring Doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

#Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = ".",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

getAvailableMatrices(proj)

proj <- filterDoublets(ArchRProj = proj)

#Dimensionality Reduction and Clustering

#ArchR implements an iterative LSI dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

#Visualizing in a 2D UMAP Embedding
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + theme(aspect.ratio=1)
png(paste0(outFile, ".cluster.umap.png"), width=3000, height=2500, res=300)
print(p1)
dev.off()

if(length(inputFiles) > 0){
  p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP") + theme(aspect.ratio=1)
  png(paste0(outFile, ".sample.umap.png"), width=3000, height=2500, res=300)
  print(p1)
  dev.off()
}

saveArchRProject(proj, outFile)

#Assigning Clusters with Gene Scores





