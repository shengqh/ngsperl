rm(list=ls()) 
sample_name='P14806_mm10'
outFile='P14806_mm10'
parSampleFile1='fileList1.txt'
parSampleFile2='fileList2.txt'
parSampleFile3=''
parFile1=''
parFile2=''
parFile3=''


setwd('/nobackup/h_cqs/shengq2/test/20260929_RNAseq_Deconvolution/BayesPrism_Deconvolution/result/P14806_mm10')

### Parameter setting end ###

library(Seurat)
#BiocManager::install("Danko-Lab/BayesPrism/BayesPrism")
library(BayesPrism)
library(logger)
library(data.table)

source("reportFunctions.R")

file_tbl=fread(parSampleFile1, data.table=FALSE, header=FALSE) |>
    dplyr::filter(V3 == sample_name)

data_options = split(file_tbl$V1, file_tbl$V2)

# change below info to your own project
ref_assay <- 'RNA'               # raw counts; never mnn.reconstructed / integrated
n_cores   <- 12                  # cores for run.prism()
# cells kept per cell type. the reference is made dense below, so memory is about
# (types x down_sample_size) cells x genes x 8 bytes - printed before it is built
down_sample_size <- 5000
pw_res    <- "."   # output directory

output_prefix = sample_name

fn_ref   = data_options$single_cell_rds  # your single cell seurat object
cell_type_column = data_options$cell_type_column    # meta.data column holding the reference cell type
fn_bulk  = data_options$count_file
gene_col = data_options$gene_column
cols_rm  = data_options$discard_columns

log_info(paste0('>>> Reading bulk RNA-seq data: ', fn_bulk))

# the delimiter is sniffed from the header line rather than taken from the extension;
# read.table's default (any whitespace) would turn a comma-separated file into one column
read_bulk <- function(fn) {
    hdr <- readLines(fn, n = 1)
    sep <- if (grepl('\t', hdr)) '\t' else if (grepl(',', hdr)) ',' else ''
    read.table(fn, header = TRUE, sep = sep, quote = '"', comment.char = '',
               stringsAsFactors = FALSE, check.names = FALSE)
}

bulk.raw <- read_bulk(fn_bulk)
if (!gene_col %in% colnames(bulk.raw)) {
    stop(fn_bulk, ": no gene column '", gene_col, "'. has: ",
            paste(utils::head(colnames(bulk.raw), 8), collapse = ", "))
}
# a column in cols_rm that is absent here is not an error, but it is reported: a
# misspelled numeric annotation column would survive and be treated as a sample
absent <- setdiff(cols_rm, colnames(bulk.raw))
if (length(absent)) {
    log_info(paste0(">>> cols_rm not present in this file (ignored): ", paste(absent, collapse = ", ")))
}
bulk.raw <- bulk.raw[, setdiff(colnames(bulk.raw), cols_rm), drop = FALSE]

# only numeric columns are samples: character annotation columns are dropped here.
# cols_rm above is for the annotation columns that ARE numeric (length, start ...)
is_num  <- vapply(bulk.raw, is.numeric, logical(1))
is_num[gene_col] <- FALSE
dropped <- setdiff(colnames(bulk.raw)[!is_num], gene_col)
if (length(dropped)) log_info(paste0(">>> dropping non-numeric columns: ", paste(dropped, collapse = ", ")))
counts <- as.matrix(bulk.raw[, is_num, drop = FALSE])
if (ncol(counts) < 2) stop("only ", ncol(counts), " numeric column(s) in ", fn_bulk)
if (any(counts < 0, na.rm = TRUE)) stop(fn_bulk, ': negative values - not a count matrix')
if (anyNA(counts)) stop(fn_bulk, ': ', sum(is.na(counts)), ' NA counts')
if (any(counts != round(counts))) log_warn(paste0(">>> WARNING: non-integer values in ", fn_bulk,
                                                     " - check it is raw counts"))
# genes without a name (failed conversion) are dropped; duplicated names are summed
genes <- trimws(as.character(bulk.raw[[gene_col]]))
keep  <- !is.na(genes) & nzchar(genes)
if (any(!keep)) log_info(paste0(">>> dropping ", sum(!keep), " rows with no gene name"))
log_info(paste0(">>> summing ", sum(duplicated(genes[keep])), " duplicated gene rows"))
bulk_matrix <- t(rowsum(counts[keep, , drop = FALSE], group = genes[keep]))
log_info(paste0(">>> Bulk matrix: ", nrow(bulk_matrix), " samples x ", ncol(bulk_matrix), " genes"))

log_info(">>> Preparing reference single-cell data ...")
# ---- reference ----
sc_object <- readRDS(fn_ref)
if (!cell_type_column %in% colnames(sc_object@meta.data)) {
    stop("No cell type column", cell_type_column, " in reference", fn_ref)
}
ref_types <- sort(unique(as.character(sc_object@meta.data[[cell_type_column]])))
log_info(paste0('--- reference cell types (', length(ref_types), ') ---'))
print(table(sc_object@meta.data[[cell_type_column]], useNA = 'ifany'))

# set idents for downsampling
Idents(sc_object) <- cell_type_column
log_info(paste0(">>> Downsampling single-cell data to ", down_sample_size, " cells per type..."))
sc_subset <- subset(sc_object, downsample = down_sample_size)
log_info(paste("Cells after downsampling:", ncol(sc_subset)))

# Seurat 5: counts split into layers (counts.1, counts.2 ...) must be joined first,
# or GetAssayData returns only one of them
if (inherits(sc_subset[[ref_assay]], 'Assay5') &&
    length(Layers(sc_subset[[ref_assay]], search = 'counts')) > 1) {
    sc_subset[[ref_assay]] <- JoinLayers(sc_subset[[ref_assay]])
}
log_info(sprintf(">>> dense reference: %d cells x %d genes, ~%.1f GB",
                 ncol(sc_subset), nrow(sc_subset),
                 ncol(sc_subset) * nrow(sc_subset) * 8 / 1e9))
sc_counts_sub <- t(as.matrix(GetAssayData(sc_subset, assay = ref_assay,
                                            slot = "counts")))
current_labels    <- as.factor(sc_subset@meta.data[[cell_type_column]])
cell_state_labels <- as.factor(sc_subset@meta.data[[cell_type_column]])
rm(sc_subset); gc()

# gene cleanup
sc.dat.filtered <- cleanup.genes(input = sc_counts_sub,
                                    input.type = "count.matrix",
                                    species = "hs",
                                    gene.group = c("Rb","Mrp","other_Rb","chrM","MALAT1"),
                                    exp.cells = 5)
rm(sc_counts_sub); gc()

# Intersect with Bulk
common_genes <- intersect(colnames(sc.dat.filtered), colnames(bulk_matrix))
log_info(paste("Common genes found:", length(common_genes)))
if (length(common_genes) < 100) stop("Too few common genes! Check formatting.")

sc.final   <- sc.dat.filtered[, common_genes]
bulk.final <- bulk_matrix[, common_genes]

log_info(">>> Creating Prism Object...")
my.prism <- new.prism(
    reference = sc.final,
    mixture = bulk.final,
    input.type = "count.matrix",
    cell.type.labels = current_labels,
    cell.state.labels = cell_state_labels,   # = the type unless a state column is given
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
)

log_info(paste0(">>> Running Gibbs Sampler (Deconvolution) on ", n_cores, " cores..."))
bp.res <- run.prism(prism = my.prism, n.cores = n_cores)

fn_rds <- paste0(output_prefix, "_BayesPrism_Object.rds")
log_info(paste0(">>> Saving BayesPrism Object: ", fn_rds))
saveRDS(bp.res, file = fn_rds)

log_info(">>> Extracting Cell Type Fractions...")
theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")
# sample ids that start with a digit come back as X1508AC1 if make.names() ran on the
# bulk header upstream; the X is removed so the ids match the original ones
theta_df <- data.frame(RNA.ID = sub("^X(?=[0-9])", "", rownames(theta), perl = TRUE),
                        as.data.frame(theta), check.names = FALSE)

rs <- rowSums(theta)
if (any(abs(rs - 1) > 0.01)) warning(output_prefix, ': fractions do not sum to 1 (',
                                        paste(round(range(rs), 3), collapse = ' .. '), ')')

fn_csv <- file.path(pw_res, paste0(output_prefix, "_Fractions.csv"))
write.csv(theta_df, fn_csv, row.names = FALSE)
log_info(paste0(">>> FINISHED: ", fn_csv, " (", nrow(theta_df), " samples x ",
        ncol(theta), " cell types)"))

