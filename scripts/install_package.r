Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")

if("BiocManager" %in% rownames(installed.packages()) == FALSE) {install.packages("BiocManager", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("hdf5r" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("hdf5r", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("remotes" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("remotes", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("BiocParallel" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("BiocParallel", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DESeq2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DESeq2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("HGNChelper" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("HGNChelper", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("MatrixEQTL" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("MatrixEQTL", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("R.utils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("R.utils", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("RColorBrewer", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Seurat" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Seurat", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SingleCellExperiment" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SingleCellExperiment", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("VennDiagram" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("VennDiagram", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("clustree" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("clustree", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("data.table" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("data.table", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("dplyr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("edgeR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("edgeR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggplot2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggrepel" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggrepel", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("grid" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("grid", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("heatmap3" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("heatmap3", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("kableExtra" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("kableExtra", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("knitr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("knitr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("lattice" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("lattice", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("mclust" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("mclust", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("readxl" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("readxl", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("reshape2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rjson" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rjson", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rmarkdown" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rmarkdown", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("scales" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scales", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("singleCellTK" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("singleCellTK", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("stringr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("stringr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tidyr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tools" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tools", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("wordcloud" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("wordcloud", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Matrix" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Matrix", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("R.utils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("R.utils", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Seurat" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Seurat", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("data.table" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("data.table", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("demuxmix" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("demuxmix", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("dnaMethyAge" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("dnaMethyAge", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("dplyr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("immunarch" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("immunarch", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("reshape2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rmarkdown" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rmarkdown", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ASCAT" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ASCAT", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ArchR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ArchR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Azimuth" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Azimuth", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("BSgenome.Hsapiens.UCSC.hg38" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("BiocParallel" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("BiocParallel", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("CHETAH" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("CHETAH", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("CellChat" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("CellChat", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ChIPpeakAnno" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ChIPpeakAnno", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ComplexHeatmap" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ComplexHeatmap", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DCATS" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DCATS", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DESeq2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DESeq2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DNAcopy" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DNAcopy", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DT" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DT", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DiffBind" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DiffBind", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DirichletReg" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DirichletReg", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DoubletFinder" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DoubletFinder", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("DropletUtils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("DropletUtils", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("EnhancedVolcano" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("EnhancedVolcano", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("EnsDb.Hsapiens.v86" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("EnsDb.Hsapiens.v86", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("GenABEL" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("GenABEL", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("GenomicRanges" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("GenomicRanges", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Hmisc" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Hmisc", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Homo.sapiens" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Homo.sapiens", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("KEGG.db" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("KEGG.db", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("KEGGprofile" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("KEGGprofile", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("MAGeCKFlute" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("MAGeCKFlute", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("cvarrichio/Matrix.utils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("cvarrichio/Matrix.utils", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("R.utils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("R.utils", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("RColorBrewer", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("RCurl" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("RCurl", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("RCy3" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("RCy3", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Rcpp" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Rcpp", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Rsamtools" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Rsamtools", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("STRINGdb" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("STRINGdb", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Seurat" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Seurat", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SeuratData" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SeuratData", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SeuratDisk" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SeuratDisk", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SeuratWrappers" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SeuratWrappers", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("Signac" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("Signac", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SignacX" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SignacX", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SingleCellExperiment" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SingleCellExperiment", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SingleR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SingleR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("SummarizedExperiment" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("SummarizedExperiment", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("TCGAbiolinks" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("TCGAbiolinks", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("TEQC" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("TEQC", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("TissueEnrich" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("TissueEnrich", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("VennDiagram" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("VennDiagram", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("WebGestaltR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("WebGestaltR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("XML" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("XML", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ape" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ape", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("arsenal" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("arsenal", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("autothresholdr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("autothresholdr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("bbmle" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("bbmle", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("biomaRt" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("biomaRt", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("biovizBase" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("biovizBase", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("celldex" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("celldex", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("choisycutoff" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("choisycutoff", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("circlize" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("circlize", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("clonevol" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("clonevol", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("cn.mops" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("cn.mops", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("colorRamps" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("colorRamps", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("cowplot" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("cowplot", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("data.table" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("data.table", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("deTS" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("deTS", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("dendextend" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("dendextend", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("dichromat" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("dichromat", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("digest" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("digest", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("dplyr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("dplyr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("edgeR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("edgeR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("flextable" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("flextable", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("forcats" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("forcats", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("formatR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("formatR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("genefilter" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("genefilter", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("getopt" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("getopt", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggExtra" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggExtra", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggbio" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggbio", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggdendro" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggdendro", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggplot2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggpubr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggpubr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggrepel" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggrepel", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggsci" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggsci", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("ggvenn" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("ggvenn", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("glmGamPoi" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("glmGamPoi", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("graph" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("graph", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("graphics" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("graphics", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("grid" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("grid", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("gridExtra", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("gtools" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("gtools", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("gtsummary" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("gtsummary", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("harmony" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("harmony", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("hdf5r" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("hdf5r", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("heatmap3" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("heatmap3", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("htmltools" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("htmltools", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("immunarch" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("immunarch", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("kableExtra" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("kableExtra", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("knitr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("knitr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("library_name, character.only = T" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("library_name, character.only = T", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("limma" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("limma", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("maftools" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("maftools", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("magrittr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("magrittr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("matrixStats" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("matrixStats", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("methylKit" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("methylKit", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("microbiome" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("microbiome", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("miloDE" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("miloDE", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("miloR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("miloR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("minfi" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("minfi", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("openxlsx" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("openxlsx", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("optparse" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("optparse", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("patchwork" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("patchwork", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("pheatmap" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("pheatmap", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("phyloseq" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("phyloseq", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("plyr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("plyr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("png" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("png", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("preprocessCore" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("preprocessCore", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("presto" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("presto", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("purrr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("purrr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("randomcoloR" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("randomcoloR", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("readr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("readr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("readxl" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("readxl", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("reshape2" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("reshape2", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rlist" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rlist", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rmarkdown" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rmarkdown", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rmdformats" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rmdformats", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rsnps" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rsnps", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rtracklayer" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rtracklayer", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("rvest" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("rvest", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("scCustomize" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scCustomize", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("scDblFinder" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scDblFinder", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("scDemultiplex" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scDemultiplex", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("scMRMA" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scMRMA", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("liuqivandy/scRNABatchQC" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("liuqivandy/scRNABatchQC", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("scales" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("scales", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("sciClone" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("sciClone", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("sparseMatrixStats" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("sparseMatrixStats", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("splicejam" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("splicejam", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("stringr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("stringr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("testit" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("testit", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("testthat" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("testthat", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tibble" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tibble", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tictoc" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tictoc", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tidyr" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tidyr", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tidyverse", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("tools" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("tools", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("umap" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("umap", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("vegan" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("vegan", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("xfun" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("xfun", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("xtable" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("xtable", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("zoo" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("zoo", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("suppressPackageStartupMessages(ComplexHeatmap" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("suppressPackageStartupMessages(ComplexHeatmap", update=FALSE, ask=FALSE, dependencies=TRUE)}
if("suppressPackageStartupMessages(ComplexHeatmap" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("suppressPackageStartupMessages(ComplexHeatmap", update=FALSE, ask=FALSE, dependencies=TRUE)}
