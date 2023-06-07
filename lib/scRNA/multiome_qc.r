rm(list=ls()) 
sample_name='AS01'
outFile='AS01'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3='fileList3.txt'
parSampleFile4='fileList4.txt'
parFile1=''
parFile2='/data/cqs/references/scrna/EnsDb.Hsapiens.v86.GRanges.UCSC.rds'
parFile3=''


setwd('/nobackup/shah_lab/shengq2/20230530_Vandy_AS_multiome/multiome_qc/result/AS01')

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

if(exists('parSampleFile5')){
  fragments_cell_map = read_file_map(parSampleFile5, do_unlist=FALSE)
}else{
  fragments_cell_map = NULL
}

counts_file = counts_file_map[[sample_name]]
if(dir.exists(counts_file)){
  counts <- Read10X(counts_file)
  atac_cells = colnames(counts)
}else{
  counts <- Read10X_h5(counts_file)
  peaks_count = counts$Peaks
  atac_cells = colnames(peaks_count)
  counts = counts$`Gene Expression`
}

# frag_cells_file = fragments_cell_map[[sample_name]]
# frag_cells = readLines(frag_cells_file)

# gex_cells = colnames(counts)
# gex_cells = gex_cells[order(gex_cells)]
# writeLines(gex_cells, paste0(outFile, ".gex.cells.txt"))
# writeLines(frag_cells, paste0(outFile, ".frag.cells.txt"))

# print(table(colnames(counts) %in% frag_cells))

# stopifnot(all(colnames(counts) %in% frag_cells))

#annotation <- readRDS(parFile2)

# create a Seurat object containing the RNA adata
obj <- CreateSeuratObject(
  counts = counts,
  assay = "RNA"
)

if(file.exists(parFile2)){
  annotation = readRDS(parFile2)
}else{
  annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotation) <- "UCSC"
  #saveRDS(annotation, 'EnsDb.Hsapiens.v86.UCSC.rds')
}

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
  cells = atac_cells
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

DefaultAssay(obj) <- "RNA"
obj <- do_sctransform(obj, use_sctransform_v2=TRUE, mc.cores=8)
obj <- RunPCA(obj)

DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

# build a joint neighbor graph using both assays
obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:20, 2:10),
  verbose = TRUE
)

obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE, graph.name="wsnn")
obj$seurat_clusters_BiMod <- obj$seurat_clusters

# build a joint UMAP visualization
obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

ctdef<-init_celltype_markers(panglao5_file = myoptions$db_markers_file,
                             species = myoptions$species,
                             curated_markers_file = myoptions$curated_markers_file,
                             HLA_panglao5_file = myoptions$HLA_panglao5_file)

cell_activity_database<-ctdef$cell_activity_database

cluster="seurat_clusters"
data_norm=get_seurat_average_expression(obj, cluster)

predict_celltype<-ORA_celltype(data_norm,cell_activity_database$cellType,cell_activity_database$weight)

if(length(predict_celltype$max_cta) > 1){
  cta_png_file=paste0(sample_name, ".cta.png")
  Plot_predictcelltype_ggplot2( predict_celltype, 
                        filename=cta_png_file)
}

new_cluster_ids<-names(predict_celltype$max_cta)
names(new_cluster_ids) = colnames(data_norm)
obj$cell_type = unlist(new_cluster_ids[as.character(obj$seurat_clusters)])
obj$seurat_cell_type = paste0(obj$seurat_clusters, ": ", obj$cell_type)

g1<-get_dim_plot_labelby(obj, reduction="umap", label.by="cell_type") + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("Cell type")
png(paste0(sample_name, ".umap.cell_type.png"), width=1800, height=1500, res=300)
print(g1)
dev.off()

g2<-get_dim_plot(obj, reduction="umap", group.by="seurat_clusters", label.by="seurat_cell_type", random_colors = FALSE) + guides(fill=guide_legend(ncol =1)) + ggtitle("Seurat cell type")
png(paste0(sample_name, ".umap.seurat_cell_type.png"), width=1800, height=1500, res=300)
print(g2)
dev.off()

if(file.exists(myoptions$bubblemap_file)){
  g<-get_bubble_plot(obj = obj, 
    bubblemap_file = myoptions$bubblemap_file, 
    group.by = "seurat_cell_type")
  dot_file = paste0(sample_name, ".dot.png")
  png(dot_file, width=get_dot_width(g), height=get_dot_height(obj, "seurat_cell_type"), res=300)
  print(g)
  dev.off()
}

saveRDS(obj, paste0(sample_name, ".obj.rds"))

# #Linking peaks to genes
# DefaultAssay(obj) <- "ATAC"

# # first compute the GC content for each peak
# obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# # link peaks to genes
# obj <- LinkPeaks(
#   object = obj,
#   peak.assay = "ATAC",
#   expression.assay = "SCT",
#   genes.use = c("COL4A5", "ADAMTS1")
# )

# p1 <- CoveragePlot(
#   object = obj,
#   region = "COL4A5",
#   features = "COL4A5",
#   expression.assay = "SCT",
#   extend.upstream = 500,
#   extend.downstream = 10000
# )

# png(paste0(sample_name, ".coverage.COL4A5.png"), width=3300, height=3000, res=300)
# print(p1)
# dev.off()
