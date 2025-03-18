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

args=commandArgs(trailingOnly=TRUE)
opt=list( input=args[1],
          output=args[2],
          title="FPKM Values Across Genomic Regions",
          width=20,
          height=10)

opt=list( input="/nobackup/brown_lab/projects/20241203_12508_cutrun_mm10/20250311_homer_overlap_diffPeaks/result/P12508_cutrun_mm10_histone.homer_overlap_diffPeaks.txt",
          output="/nobackup/brown_lab/projects/20241203_12508_cutrun_mm10/20250311_homer_overlap_diffPeaks_vis/result/P12508_cutrun_mm10_histone.homer_overlap_diffPeaks.png",
          title="FPKM Values Across Genomic Regions",
          width=20,
          height=10)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be specified", call.=FALSE)
}

# Read the data
cat("Reading data from", opt$input, "...\n")
data <- read.delim(opt$input, stringsAsFactors=FALSE, check.names=FALSE) |>
  dplyr::rename(Chr=bed_chr, Start=bed_start, End=bed_end, FPKM=sample_fpkm, Sample=sample_name)
# Check if data was loaded correctly
if (nrow(data) == 0) {
  stop("No data found in the input file")
}

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
sorted_data <- data %>%
  # Sort by chromosome (naturally) and start position
  arrange(chr_order(Chr), Start)

sorted_data$Chr=factor(sorted_data$Chr, levels=unique(sorted_data$Chr))

# Create the plot
cat("Creating plot...\n")

# Create the plot
p <- ggplot(sorted_data, aes(x = Start, y = FPKM)) +
  geom_point() +
  theme_bw() +
  labs(title = opt$title,
       x = "Genomic Position",
       y = "FPKM Value") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_blank()) +
  facet_grid(Sample ~ Chr, scales = "free_x", space = "free_x")

# Save the plot
cat("Saving plot to", opt$output, "...\n")
ggsave(opt$output, plot = p, width = 12, height = 7, units="in", dpi=300, bg="white")

cat("Done!\n")


