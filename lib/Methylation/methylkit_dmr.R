suppressPackageStartupMessages({
  library(Cairo)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(ggplot2)
  library(methylKit)
  library(rtracklayer)
})

read_parameter_map <- function(filename) {
  values <- read.table(filename, sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (ncol(values) < 2) {
    stop("Parameter file must contain value and name columns: ", filename)
  }
  result <- as.list(values[[1]])
  names(result) <- values[[2]]
  result
}

get_parameter <- function(params, name, default = NULL, numeric = FALSE) {
  value <- params[[name]]
  if (is.null(value) || is.na(value) || value == "") {
    value <- default
  }
  if (numeric && !is.null(value)) {
    value <- as.numeric(value)
  }
  value
}

require_nonempty_file <- function(filename, description) {
  if (is.null(filename) || filename == "" || !file.exists(filename) || file.info(filename)$size == 0) {
    stop(description, " is missing or empty: ", filename)
  }
}

assign_gene_context <- function(target, gene_features) {
  context <- rep("intergenic", length(target))
  feature_order <- c(introns = "intron", exons = "exon", promoters = "promoter")
  for (feature_name in names(feature_order)) {
    hits <- findOverlaps(target, gene_features[[feature_name]], ignore.strand = TRUE)
    context[unique(queryHits(hits))] <- feature_order[[feature_name]]
  }
  context
}

assign_cpg_context <- function(target, islands, shores) {
  context <- rep("other", length(target))
  shore_hits <- findOverlaps(target, shores, ignore.strand = TRUE)
  context[unique(queryHits(shore_hits))] <- "shore"
  island_hits <- findOverlaps(target, islands, ignore.strand = TRUE)
  context[unique(queryHits(island_hits))] <- "CpG_island"
  context
}

read_bed12_gene_features <- function(filename, promoter_up, promoter_down) {
  bed <- data.table::fread(filename, header = FALSE, sep = "\t", fill = TRUE,
                          data.table = FALSE, showProgress = FALSE)
  if (ncol(bed) < 12) {
    stop("Transcript BED annotation must contain at least 12 columns: ", filename)
  }

  transcript_names <- as.character(bed[[4]])
  transcript_names[is.na(transcript_names) | transcript_names == ""] <-
    paste0("transcript_", which(is.na(transcript_names) | transcript_names == ""))
  strands <- as.character(bed[[6]])
  if (any(!strands %in% c("+", "-"))) stop("Transcript BED contains invalid or missing strands")
  tx <- GRanges(seqnames = bed[[1]],
                ranges = IRanges(start = as.integer(bed[[2]]) + 1L, end = as.integer(bed[[3]])),
                strand = strands, transcript = transcript_names)

  exon_list <- vector("list", nrow(bed))
  intron_list <- vector("list", nrow(bed))
  for (i in seq_len(nrow(bed))) {
    block_count <- as.integer(bed[[10]][i])
    block_sizes <- as.integer(Filter(nzchar, strsplit(as.character(bed[[11]][i]), ",", fixed = TRUE)[[1]]))
    block_starts <- as.integer(Filter(nzchar, strsplit(as.character(bed[[12]][i]), ",", fixed = TRUE)[[1]]))
    if (is.na(block_count) || block_count < 1 || length(block_sizes) != block_count ||
        length(block_starts) != block_count) {
      stop("Invalid BED12 exon blocks for transcript ", transcript_names[i])
    }
    exon_starts <- as.integer(bed[[2]][i]) + block_starts + 1L
    exon_ends <- exon_starts + block_sizes - 1L
    exon_list[[i]] <- GRanges(seqnames = bed[[1]][i],
                              ranges = IRanges(exon_starts, exon_ends),
                              strand = strands[i], transcript = transcript_names[i])
    if (block_count > 1) {
      intron_starts <- exon_ends[-block_count] + 1L
      intron_ends <- exon_starts[-1] - 1L
      keep <- intron_starts <= intron_ends
      intron_list[[i]] <- GRanges(seqnames = bed[[1]][i],
                                  ranges = IRanges(intron_starts[keep], intron_ends[keep]),
                                  strand = strands[i], transcript = transcript_names[i])
    } else {
      intron_list[[i]] <- GRanges()
    }
  }

  tss <- resize(tx, width = 1, fix = "start")
  promoters <- GenomicRanges::promoters(tx, upstream = as.integer(promoter_up),
                                        downstream = as.integer(promoter_down))
  list(exons = unlist(GRangesList(exon_list), use.names = FALSE),
       introns = unlist(GRangesList(intron_list), use.names = FALSE),
       promoters = promoters, tss = tss)
}

read_gene_features <- function(filename, promoter_up, promoter_down) {
  extension <- tolower(sub("[.]gz$", "", filename))
  is_gff <- grepl("[.](gtf|gff|gff3)$", extension)
  if (!is_gff) return(read_bed12_gene_features(filename, promoter_up, promoter_down))

  txdb <- GenomicFeatures::makeTxDbFromGFF(filename)
  tx <- GenomicFeatures::transcripts(txdb, columns = c("tx_id", "tx_name"))
  tx_name <- as.character(mcols(tx)$tx_name)
  tx_name[is.na(tx_name) | tx_name == ""] <- as.character(mcols(tx)$tx_id[is.na(tx_name) | tx_name == ""])
  mcols(tx)$transcript <- tx_name
  exon_ranges <- GenomicFeatures::exons(txdb, columns = c("tx_name"))
  intron_ranges <- unlist(GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE), use.names = TRUE)
  promoters <- GenomicRanges::promoters(tx, upstream = as.integer(promoter_up),
                                        downstream = as.integer(promoter_down))
  tss <- resize(tx, width = 1, fix = "start")
  list(exons = exon_ranges, introns = intron_ranges, promoters = promoters, tss = tss)
}

read_cpg_features <- function(filename, shore_size) {
  islands <- rtracklayer::import(filename, format = "BED")
  islands <- reduce(islands, ignore.strand = TRUE)
  strand(islands) <- "*"
  upstream <- flank(islands, width = as.integer(shore_size), start = TRUE, both = FALSE)
  downstream <- flank(islands, width = as.integer(shore_size), start = FALSE, both = FALSE)
  shores <- setdiff(reduce(c(upstream, downstream), ignore.strand = TRUE), islands,
                    ignore.strand = TRUE)
  list(islands = islands, shores = shores)
}

annotation_summary <- function(values, annotation_name, levels) {
  counts <- table(factor(values, levels = levels))
  total <- sum(counts)
  data.frame(
    annotation = annotation_name,
    category = names(counts),
    count = as.integer(counts),
    percent = if (total == 0) rep(0, length(counts)) else as.numeric(counts) * 100 / total,
    stringsAsFactors = FALSE
  )
}

params <- read_parameter_map(parSampleFile2)
assembly <- get_parameter(params, "assembly")
pipeline <- get_parameter(params, "pipeline", "amp")
mincov <- get_parameter(params, "mincov", 3, numeric = TRUE)
high_cov_pct <- get_parameter(params, "high_cov_pct", 99.99, numeric = TRUE)
window_size <- get_parameter(params, "window_size", 1000, numeric = TRUE)
step_size <- get_parameter(params, "step_size", 1000, numeric = TRUE)
min_cpgs <- get_parameter(params, "min_cpgs", 10, numeric = TRUE)
min_per_group <- get_parameter(params, "min_per_group", 0, numeric = TRUE)
difference <- get_parameter(params, "difference", 25, numeric = TRUE)
qvalue <- get_parameter(params, "qvalue", 0.01, numeric = TRUE)
ncore <- get_parameter(params, "ncore", 1, numeric = TRUE)
overdispersion <- get_parameter(params, "overdispersion", "MN")
test_method_key <- tolower(get_parameter(params, "test_method", "dss"))
test_method <- switch(test_method_key,
                      dss = "dss",
                      f = "F",
                      chisq = "Chisq",
                      fast.fisher = "fast.fisher",
                      midpval = "midPval",
                      stop("Unsupported methylKit DMR test_method: ", test_method_key))
adjust <- get_parameter(params, "adjust", "BH")
use_raw_pvalue <- get_parameter(params, "use_raw_pvalue", 0, numeric = TRUE)
promoter_up <- get_parameter(params, "promoter_up", 1000, numeric = TRUE)
promoter_down <- get_parameter(params, "promoter_down", 1000, numeric = TRUE)
shore_size <- get_parameter(params, "shore_size", 2000, numeric = TRUE)

if (is.null(assembly) || assembly == "") stop("assembly/genome is not defined")
if (window_size < 1 || step_size < 1 || min_cpgs < 1) stop("Window, step, and minimum CpG values must be positive")
if (mincov < 0 || high_cov_pct < 0 || high_cov_pct > 100 ||
    promoter_up < 0 || promoter_down < 1 || shore_size < 1) {
  stop("Coverage and annotation flank sizes are invalid")
}
comparisons <- read.table(parSampleFile3, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
comparison <- comparisons[comparisons[[2]] == sample_name, , drop = FALSE]
if (nrow(comparison) != 2) {
  stop("Comparison ", sample_name, " must contain exactly two group names; found ", nrow(comparison))
}
control_group_name <- comparison[[1]][1]
treatment_group_name <- comparison[[1]][2]

groups <- read.table(parSampleFile4, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(groups) < 2) stop("Group definition must contain sample and group columns")
control_names <- unique(groups[[1]][groups[[2]] == control_group_name])
treatment_names <- unique(groups[[1]][groups[[2]] == treatment_group_name])
if (length(control_names) == 0 || length(treatment_names) == 0) {
  stop("Both groups must contain samples. Control=", control_group_name,
       "; treatment=", treatment_group_name)
}
if (min_per_group > min(length(control_names), length(treatment_names))) {
  stop("min_per_group cannot exceed the number of samples in the smaller group")
}
samples <- c(control_names, treatment_names)
treatment <- c(rep(0, length(control_names)), rep(1, length(treatment_names)))

input_map <- read.table(parSampleFile5, sep = "\t", header = FALSE,
                        stringsAsFactors = FALSE, quote = "", comment.char = "")
if (ncol(input_map) < 2) stop("MethylKit input list must contain file and sample columns")
colnames(input_map)[1:2] <- c("file", "sample")
duplicated_samples <- unique(input_map$sample[duplicated(input_map$sample)])
if (length(duplicated_samples) > 0) {
  stop("Multiple prepared CpG files found for sample(s): ", paste(duplicated_samples, collapse = ", "))
}
missing_samples <- setdiff(samples, input_map$sample)
if (length(missing_samples) > 0) {
  stop("No prepared CpG input found for sample(s): ", paste(missing_samples, collapse = ", "))
}
input_map <- input_map[match(samples, input_map$sample), , drop = FALSE]
bad_files <- input_map$file[!file.exists(input_map$file) | file.info(input_map$file)$size == 0]
if (length(bad_files) > 0) {
  stop("Prepared CpG input file(s) missing or empty: ", paste(bad_files, collapse = ", "))
}
require_nonempty_file(parFile1, "Transcript BED annotation")
require_nonempty_file(parFile2, "CpG island BED annotation")

cat("Comparison:", treatment_group_name, "vs", control_group_name, "\n")
cat("Samples:", paste(samples, collapse = ", "), "\n")
cat("Tiling windows:", window_size, "bp; step:", step_size,
    "bp; minimum CpGs:", min_cpgs, "\n")

cpg_obj <- methRead(
  location = as.list(input_map$file),
  sample.id = as.list(samples),
  assembly = assembly,
  treatment = treatment,
  context = "CpG",
  mincov = mincov,
  pipeline = pipeline
)
if (high_cov_pct > 0 && high_cov_pct < 100) {
  cpg_obj <- filterByCoverage(cpg_obj, lo.count = mincov, lo.perc = NULL,
                              hi.count = NULL, hi.perc = high_cov_pct)
} else {
  cpg_obj <- filterByCoverage(cpg_obj, lo.count = mincov, lo.perc = NULL,
                              hi.count = NULL, hi.perc = NULL)
}
tiles <- tileMethylCounts(cpg_obj, win.size = as.integer(window_size),
                          step.size = as.integer(step_size), cov.bases = as.integer(min_cpgs))
rm(cpg_obj)

if (min_per_group > 0) {
  tiled_meth <- unite(tiles, destrand = FALSE, min.per.group = as.integer(min_per_group))
} else {
  tiled_meth <- unite(tiles, destrand = FALSE)
}
rm(tiles)
if (nrow(getData(tiled_meth)) == 0) stop("No tiled regions passed coverage requirements for all required samples")

if (test_method == "dss") {
  diff_all <- calculateDiffMethDSS(tiled_meth, adjust = adjust, mc.cores = ncore)
} else {
  diff_all <- calculateDiffMeth(tiled_meth, overdispersion = overdispersion,
                                adjust = adjust, test = test_method, mc.cores = ncore)
}
rm(tiled_meth)

prefix <- paste0(sample_name, ".methylkit.dmr")
saveRDS(diff_all, paste0(prefix, ".rds"))

if (use_raw_pvalue != 0) {
  diff_all$adjusted_qvalue <- diff_all$qvalue
  diff_all$qvalue <- diff_all$pvalue
}
dmrs <- getMethylDiff(diff_all, difference = difference, qvalue = qvalue, type = "all")
dmr_df <- as.data.frame(dmrs)
if (use_raw_pvalue != 0 && nrow(dmr_df) > 0) {
  dmr_df$qvalue <- dmr_df$adjusted_qvalue
  dmr_df$adjusted_qvalue <- NULL
}
dmr_df$direction <- ifelse(dmr_df$meth.diff > 0,
                           paste0("hyper_in_", treatment_group_name),
                           paste0("hyper_in_", control_group_name))
if (nrow(dmr_df) > 0) dmr_df <- dmr_df[order(dmr_df$qvalue), , drop = FALSE]
write.table(dmr_df, paste0(sample_name, ".methylkit.dmrs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

sample_table <- data.frame(
  sample = samples,
  group = c(rep(control_group_name, length(control_names)),
            rep(treatment_group_name, length(treatment_names))),
  role = c(rep("control", length(control_names)), rep("treatment", length(treatment_names))),
  input_file = input_map$file,
  stringsAsFactors = FALSE
)
write.table(sample_table, paste0(prefix, ".samples.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

parameter_table <- data.frame(
  parameter = c("comparison", "control_group", "treatment_group", "assembly", "window_size",
                "step_size", "min_cpgs", "mincov", "high_cov_pct", "min_per_group",
                "test_method", "overdispersion", "adjust", "difference", "qvalue",
                "use_raw_pvalue", "promoter_up", "promoter_down", "shore_size",
                "tested_regions", "significant_dmrs"),
  value = c(sample_name, control_group_name, treatment_group_name, assembly, window_size,
            step_size, min_cpgs, mincov, high_cov_pct, min_per_group, test_method,
            overdispersion, adjust, difference, qvalue, use_raw_pvalue, promoter_up,
            promoter_down, shore_size,
            nrow(getData(diff_all)), nrow(dmr_df)),
  stringsAsFactors = FALSE
)
write.table(parameter_table, paste0(prefix, ".parameters.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

gene_levels <- c("promoter", "exon", "intron", "intergenic")
cpg_levels <- c("CpG_island", "shore", "other")
if (nrow(dmr_df) > 0) {
  dmr_ranges <- as(dmrs, "GRanges")
  gene_features <- read_gene_features(parFile1, promoter_up, promoter_down)
  dmr_df$gene_context <- assign_gene_context(dmr_ranges, gene_features)

  dmr_df$nearest_transcript <- NA_character_
  dmr_df$distance_to_tss <- NA_real_
  nearest_tss <- distanceToNearest(dmr_ranges, gene_features$tss, ignore.strand = TRUE)
  if (length(nearest_tss) > 0) {
    target_rows <- queryHits(nearest_tss)
    tss_rows <- subjectHits(nearest_tss)
    dmr_df$nearest_transcript[target_rows] <- as.character(mcols(gene_features$tss)$transcript[tss_rows])
    dmr_df$distance_to_tss[target_rows] <- mcols(nearest_tss)$distance
  }

  cpg_features <- read_cpg_features(parFile2, shore_size)
  dmr_df$cpg_context <- assign_cpg_context(dmr_ranges, cpg_features$islands, cpg_features$shores)
} else {
  dmr_df$gene_context <- character(0)
  dmr_df$nearest_transcript <- character(0)
  dmr_df$distance_to_tss <- numeric(0)
  dmr_df$cpg_context <- character(0)
}
write.table(dmr_df, paste0(prefix, ".annotated.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

summary_table <- rbind(
  annotation_summary(dmr_df$gene_context, "gene_context", gene_levels),
  annotation_summary(dmr_df$cpg_context, "cpg_context", cpg_levels)
)
write.table(summary_table, paste0(prefix, ".annotation_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

if (nrow(dmr_df) > 0) {
  chromosome_summary <- as.data.frame(table(chromosome = dmr_df$chr, direction = dmr_df$direction),
                                      stringsAsFactors = FALSE)
  chromosome_summary <- chromosome_summary[chromosome_summary$Freq > 0, , drop = FALSE]
  colnames(chromosome_summary)[3] <- "count"
} else {
  chromosome_summary <- data.frame(chromosome = character(), direction = character(), count = integer())
}
write.table(chromosome_summary, paste0(prefix, ".chromosome_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

plot_annotation <- function(summary_data, annotation_name, output_file, title) {
  plot_data <- summary_data[summary_data$annotation == annotation_name, , drop = FALSE]
  graph <- ggplot(plot_data, aes(x = category, y = percent, fill = category)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = count), vjust = -0.25, size = 3.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "DMRs (%)", title = title) +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave(output_file, graph, width = 6, height = 4, dpi = 300, bg = "white")
}

plot_annotation(summary_table, "gene_context", paste0(prefix, ".annotation.png"),
                "DMR gene-context annotation")
plot_annotation(summary_table, "cpg_context", paste0(prefix, ".cpg_context.png"),
                "DMR CpG-island context")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
