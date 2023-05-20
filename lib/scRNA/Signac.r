rm(list=ls()) 
sample_name='AS01'
outFile='AS01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parFile1='/data/cqs/references/scrna/EnsDb.Hsapiens.v86.GRanges.UCSC.rds'
parFile2=''
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20230517_Vandy_AS_scRNA_sct2_atac/signac_individual/result/AS01')

### Parameter setting end ###

source("scRNA_func.r")
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(data.table)
library(biovizBase)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

options(future.globals.maxSize= 10779361280)
random.seed=20200107

myoptions<-read_file_map(parSampleFile2, do_unlist=FALSE)
nCount_cutoff=to_numeric(myoptions$nCount_cutoff, 500)
nCount_cutoff_max=to_numeric(myoptions$nCount_cutoff_max, 200000)
nFeature_cutoff_min=to_numeric(myoptions$nFeature_cutoff_min, 200)
nFeature_cutoff_max=to_numeric(myoptions$nFeature_cutoff_max, 100000)
nCount_ATAC_min=to_numeric(myoptions$nCount_ATAC_min, 1000)
nCount_ATAC_max=to_numeric(myoptions$nCount_ATAC_max, 100000)
max_nucleosome_signal=to_numeric(myoptions$max_nucleosome_signal, 2)
min_TSS_enrichment=to_numeric(myoptions$min_TSS_enrichment, 1)
macs_path=get_value(myoptions$macs_path, NULL)

counts_file_map = read_file_map(parSampleFile1, do_unlist=FALSE)
fragments_file_map = read_file_map(parSampleFile3, do_unlist=FALSE)

counts_file = counts_file_map[[sample_name]]
counts <- Read10X_h5(counts_file)
fragpath = fragments_file_map[[sample_name]]

annotation <- readRDS(parFile1)

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

png(paste0(outFile, ".atac.qc.png"), width=4000, height=1000, res=300)
g<-VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
) * xlab("") * theme(axis.ticks.x=element_blank(),
                     axis.text.x = element_blank())
print(g)
dev.off()

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < nCount_ATAC_max &
    nCount_RNA < nCount_cutoff_max &
    nCount_ATAC > nCount_ATAC_min &
    nCount_RNA > nCount_cutoff &
    nucleosome_signal < max_nucleosome_signal &
    TSS.enrichment > min_TSS_enrichment
)

saveRDS(pbmc, "obj.rds")

frags <- Fragments(object = pbmc)
allfragpaths <- sapply(X = frags, FUN = GetFragmentData, slot = "path")
fragpath <- GetFragmentData(object = allfragpaths, slot = "path")

peaks <- CallPeaks(pbmc, macs2.path=macs_path)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
