#!/usr/bin/env Rscript

# Script to visualize FPKM values across genomic regions for multiple samples
# Input: A tab-delimited file with columns: chr, start, end, SAMPLE1_FPKM, SAMPLE2_FPKM, ..., SAMPLEN_FPKM
# Output: PDF file with line plot showing FPKM values across genomic regions for each sample

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(scales)  # For pretty_breaks
})

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input tab-delimited file", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="fpkm_plot.pdf", 
              help="Output PDF file [default= %default]", metavar="FILE"),
  make_option(c("-t", "--title"), type="character", default="FPKM Values Across Genomic Regions", 
              help="Plot title [default= %default]", metavar="STRING"),
  make_option(c("-w", "--width"), type="numeric", default=20, 
              help="Plot width in inches [default= %default]", metavar="NUMBER"),
  make_option(c("-h", "--height"), type="numeric", default=10, 
              help="Plot height in inches [default= %default]", metavar="NUMBER")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

opt=list( input="/nobackup/brown_lab/projects/20241203_12508_cutrun_mm10/overlap/result/P12508_cutrun_mm10_histone.overlap.txt",
          output="/nobackup/brown_lab/projects/20241203_12508_cutrun_mm10/overlap/result/P12508_cutrun_mm10_histone.overlap.txt.png",
          title="FPKM Values Across Genomic Regions",
          width=20,
          height=10)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified", call.=FALSE)
}

# Read the data
cat("Reading data from", opt$input, "...\n")
data <- read.delim(opt$input, stringsAsFactors=FALSE, check.names=FALSE)

# Check if data was loaded correctly
if (nrow(data) == 0) {
  stop("No data found in the input file")
}

# Column names should be Chr, Start, End followed by sample names
if (!all(c("Chr", "Start", "End") %in% colnames(data))) {
  stop("Input file must have columns named 'Chr', 'Start', and 'End'")
}

# Get sample columns (all columns except Chr, Start, End)
sample_cols <- setdiff(colnames(data), c("Chr", "Start", "End"))
if (length(sample_cols) == 0) {
  stop("No sample columns found in the input file")
}

cat("Found", length(sample_cols), "samples:", paste(sample_cols, collapse=", "), "\n")

# Function to standardize chromosome sorting
chr_order <- function(chr_vec) {
  # Extract numeric part from chromosome names
  num_part <- gsub("chr", "", chr_vec)
  # Replace X, Y with high numbers to sort after numeric chromosomes
  num_part <- gsub("X", "23", num_part)
  num_part <- gsub("Y", "24", num_part)
  num_part <- gsub("M", "25", num_part)
  # Convert to numeric for proper sorting
  as.numeric(num_part)
}

# Convert to long format and sort by chromosome and start position
data_long <- data %>%
  # Sort by chromosome (naturally) and start position
  arrange(chr_order(Chr), Start) %>%
  # Convert to long format
  pivot_longer(cols = all_of(sample_cols),
               names_to = "Sample", 
               values_to = "FPKM") %>%
  # Convert FPKM to numeric
  mutate(FPKM = as.numeric(FPKM))

# Create a continuous genomic coordinate system
# First, sort data by chromosome and position
sorted_data <- data_long %>%
  arrange(chr_order(Chr), Start)

# Create chromosome boundaries for plotting
chrom_data <- sorted_data %>%
  group_by(Chr) %>%
  summarise(start_pos = min(Start),
            end_pos = max(End),
            .groups = 'drop')

# Add cumulative position to create continuous x-axis
chrom_data <- chrom_data %>%
  mutate(length = end_pos - start_pos) %>%
  arrange(chr_order(Chr)) %>%
  mutate(cum_start = cumsum(c(0, head(length, -1))),
         cum_end = cumsum(length))

# Join to get cumulative position for each region
data_long <- data_long %>%
  left_join(chrom_data %>% select(Chr, cum_start), by = "Chr") %>%
  mutate(abs_pos = Start + cum_start)

# for fpkm, replace all extremely higher values to 99% quantile
data_long$FPKM <- ifelse(data_long$FPKM > quantile(data_long$FPKM, 0.999), quantile(data_long$FPKM, 0.999), data_long$FPKM)

# Create the plot
cat("Creating plot...\n")
p <- ggplot(data_long, aes(x = abs_pos, y = FPKM, group = interaction(Sample, Chr), color = Sample)) +
  geom_point(size = 1) +
  theme_bw() +
  labs(title = opt$title,
       x = "Genomic Position",
       y = "FPKM Value") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(
    breaks = chrom_data$cum_start + (chrom_data$length / 2),
    labels = chrom_data$Chr
  ) +
  theme(panel.grid.minor.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

# Add chromosome boundaries
p <- p + geom_vline(xintercept = chrom_data$cum_end[-nrow(chrom_data)], 
                    linetype = "dashed", 
                    color = "gray70", 
                    alpha = 0.5) + facet_grid(Sample~.)

# Save the plot
cat("Saving plot to", opt$output, "...\n")
ggsave(opt$output, plot = p, width = opt$width, height = opt$height, units="in", dpi=300, bg="white")

cat("Done!\n")


