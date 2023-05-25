rm(list=ls()) 
outFile='as_multiome'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20230517_Vandy_AS_scRNA_sct2_atac/archR_all/result')

### Parameter setting end ###

source("scRNA_func.r")
#https://www.archrproject.com/articles/Articles/tutorial.html
#https://github.com/GreenleafLab/ArchR_Website/blob/ui_updates/bookdown/21_MultiomicDataAnalysis.Rmd

#https://www.archrproject.com/bookdown/index.html#section

library(ArchR)
library(presto)

set.seed(20230512)

myoptions<-read_file_map("fileList2.txt", do_unlist=FALSE)
minTSS=to_numeric(myoptions$minTSS, 4)
minFrags=to_numeric(myoptions$minFrags, 1000)
threads=to_numeric(myoptions$threads, 8)
genome=get_value(myoptions$genome, "hg38")

lsi_resolution=to_numeric(myoptions$lsi_resolution, 0.2)
lsi_sampleCells=to_numeric(myoptions$lsi_sampleCells, 10000)

cluster_resolution=to_numeric(myoptions$cluster_resolution, 0.8)

umap_minDist=to_numeric(myoptions$umap_minDist, 0.5)

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

getTSSEnrichmentPlot<-function(Metadata, sample_name, filterTSS, filterFrags) {
  ggtitle <- sprintf("%s\n%s\n%s",
                      paste0(sample_name, "\nnCells Pass Filter = ", sum(Metadata$Keep)),
                      paste0("Median Frags = ", median(Metadata$nFrags[Metadata$Keep==1])),
                      paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment[Metadata$Keep==1]))
  )
  
  gg <- ggPoint(
    x = pmin(log10(Metadata$nFrags), 5) + rnorm(length(Metadata$nFrags), sd = 0.00001),
    y = Metadata$TSSEnrichment + rnorm(length(Metadata$nFrags), sd = 0.00001),
    colorDensity = TRUE,
    xlim = c(2.5, 5),
    ylim = c(0, max(Metadata$TSSEnrichment) * 1.05),
    baseSize = 6,
    continuousSet = "sambaNight",
    xlabel = "Log 10 (Unique Fragments)",
    ylabel = "TSS Enrichment",
    title = ggtitle,
    rastr = TRUE) +
    geom_hline(yintercept=filterTSS, lty = "dashed", linewidth = 0.25) +
    geom_vline(xintercept=log10(filterFrags), lty = "dashed", linewidth = 0.25)
  
  return(gg)
}

g = NULL
for(s_name in names(inputFiles)){
  Metadata<-readRDS(paste0('QualityControl/', s_name, "/", s_name, "-Pre-Filter-Metadata.rds"))
  gtss<-getTSSEnrichmentPlot(Metadata, s_name, minTSS, minFrags) + theme(aspect.ratio=1)
  gtss$data$Sample=s_name
  if(is.null(g)){
    g=gtss
  }else{
    g$data = rbind(g$data, gtss$data)
  }
}
nsample=length(inputFiles)
ncol=ceiling(sqrt(nsample))
nrow=ceiling(nsample/ncol)
gg<-g+facet_wrap(~Sample) + theme_bw3() + ggtitle("") + NoLegend()
png(paste0(outFile, ".TSSEnrichment.png"), width=ncol * 600, height=nrow*600, res=300)
print(gg)
dev.off()

#Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = ".",
  copyArrows = FALSE #
  #copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#Filtering doublets
proj <- addDoubletScores(
  input = proj, 
  useMatrix = "TileMatrix", 
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1)

proj <- filterDoublets(proj)

getAvailableMatrices(proj)

#ArchR implements an iterative LSI dimensionality reduction
proj <- addIterativeLSI(
  ArchRProj = proj, 
  useMatrix = "TileMatrix", 
  name = "LSI_ATAC",
  iterations = 2, 
  clusterParams = list(
    resolution = lsi_resolution, 
    sampleCells = lsi_sampleCells,
    n.start = 10
  ),
  saveIterations = FALSE,
  varFeatures = 25000, 
  dimsToUse = 1:30  
)

#6.1 Uniform Manifold Approximation and Projection (UMAP)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "LSI_ATAC", 
  name = "UMAP_ATAC", 
  nNeighbors = 30, 
  minDist = umap_minDist, 
  metric = "cosine"
)

#5.1 Clustering using Seurat’s FindClusters() function
proj <- addClusters(
  input = proj,
  reducedDims = "LSI_ATAC",
  method = "Seurat",
  name = "Clusters_ATAC",
  resolution = cluster_resolution
)

table(proj$Clusters_ATAC)

p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "cellColData", 
  name = "Clusters_ATAC", 
  embedding = "UMAP_ATAC") + theme(aspect.ratio=1)
png(paste0(outFile, ".cluster.umap.png"), width=2000, height=1700, res=300)
print(p1)
dev.off()

p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "cellColData", 
  name = "Sample", 
  embedding = "UMAP_ATAC") + theme(aspect.ratio=1) 
p1 <- p1 + facet_wrap(~color) + theme_bw3() + scale_color_manual(values=rep("black", nsample)) + NoLegend()
png(paste0(outFile, ".sample.umap.all.png"), width=ncol * 600, height=nrow*600, res=300)
print(p1)
dev.off()

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters_ATAC")
proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "Clusters_ATAC")
proj <- addPeakMatrix(ArchRProj = proj)

rds_file = paste0(outFile, ".atac_peaks.rds")
saveRDS(proj, rds_file)

# 
# #7.1 Calculating Gene Scores in ArchR
# proj <- addGeneScoreMatrix(proj, force=TRUE)
# 
# #7.3 Identifying Marker Genes
# markersGS <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = "GeneScoreMatrix", 
#   groupBy = "Clusters_ATAC",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# cutOff="FDR <= 0.01 & Log2FC >= 1.25"
# 
# markerList <- getMarkers(
#   markersGS, 
#   cutOff = cutOff)
# 
# g <- plotMarkerHeatmap(
#   seMarker = markersGS, 
#   cutOff = cutOff, 
#   transpose = TRUE
# )
# png(paste0(outFile, ".markers.heatmap.png"), width=3000, height=1000, res=300)
# print(g)
# dev.off()
# 
# #7.5 Marker Genes Imputation with MAGIC, for visualization
# proj <- addImputeWeights(proj, reducedDims = "LSI_ATAC")
# 

