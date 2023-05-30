rm(list=ls()) 
sample_name='ENCDO068KYD'
outFile='ENCDO068KYD'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parSampleFile5='fileList5.txt'
parFile1=''
parFile2='/data/cqs/references/scrna/EnsDb.Hsapiens.v86.GRanges.UCSC.rds'
parFile3=''


setwd('/nobackup/brown_lab/projects/20230523_encode_liver_multiome/multiome_qc/result/ENCDO068KYD')

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
macs2_path=get_value(myoptions$macs2_path, NULL)

counts_file_map = read_file_map(parSampleFile1, do_unlist=FALSE)
fragments_file_map = read_file_map(parSampleFile3, do_unlist=FALSE)
fragments_cell_map = read_file_map(parSampleFile5, do_unlist=FALSE)

counts_file = counts_file_map[[sample_name]]
if(dir.exists(counts_file)){
  counts <- Read10X(counts_file)
}else{
  counts <- Read10X_h5(counts_file)
}

frag_cells_file = fragments_cell_map[[sample_name]]
frag_cells = readLines(frag_cells_file)

gex_cells = colnames(counts)
gex_cells = gex_cells[order(gex_cells)]
writeLines(gex_cells, paste0(outFile, ".gex.cells.txt"))
writeLines(frag_cells, paste0(outFile, ".frag.cells.txt"))

print(table(colnames(counts) %in% frag_cells))

stopifnot(all(colnames(counts) %in% frag_cells))

#annotation <- readRDS(parFile2)

# create a Seurat object containing the RNA adata
obj <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
saveRDS(annotation, 'EnsDb.Hsapiens.v86.UCSC.rds')

fragpath = fragments_file_map[[sample_name]]

# create fragment object
frag_obj = CreateFragmentObject(fragpath)



peaks = CallPeaks(frag_obj)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = frag_obj,
  features = peaks,
  cells = colnames(obj)
)

# create ATAC assay and add it to the object
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(obj) <- "ATAC"

obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

png(paste0(outFile, ".atac.qc.png"), width=4000, height=1000, res=300)
g<-VlnPlot(
  object = obj,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
) * xlab("") * theme(axis.ticks.x=element_blank(),
                     axis.text.x = element_blank())
print(g)
dev.off()

# filter out low quality cells
obj <- subset(
  x = obj,
  subset = nCount_ATAC < nCount_ATAC_max &
    nCount_RNA < nCount_cutoff_max &
    nCount_ATAC > nCount_ATAC_min &
    nCount_RNA > nCount_cutoff &
    nucleosome_signal < max_nucleosome_signal &
    TSS.enrichment > min_TSS_enrichment
)

saveRDS(obj, "obj.rds")

frags <- Fragments(object = obj)
allfragpaths <- sapply(X = frags, FUN = GetFragmentData, slot = "path")
fragpath <- GetFragmentData(object = allfragpaths, slot = "path")

peaks <- CallPeaks(obj, macs2.path=macs_path)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(obj),
  features = peaks,
  cells = colnames(obj)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
