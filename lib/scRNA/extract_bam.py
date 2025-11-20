import pysam
import argparse
import sys
import logging
import os

DEBUG = False
NOT_DEBUG = not DEBUG

parser = argparse.ArgumentParser(description="Extract BAM reads for specific cells.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input_bam', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
parser.add_argument('-c', '--cell_file', action='store', nargs='?', help='Input cell file, one cell barcode per line', required=NOT_DEBUG)
parser.add_argument('-o', '--output_bam', action='store', nargs='?', help="Output BAM file", required=NOT_DEBUG)
parser.add_argument('-b', '--cell_barcode_tag', action='store', nargs='?', default='CB', help='Cell barcode tag in BAM file')
parser.add_argument('--index', action='store_true', help='Index the output BAM file after creation')
args = parser.parse_args()

if DEBUG:
  args.input_bam = "/path/to/input.bam"
  args.cell_file = "/path/to/cells.txt"
  args.output_bam = "/path/to/output.bam"
  args.cell_barcode_tag = "CB"
  args.index = True

logger = logging.getLogger('extract_bam')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

# Read cell barcodes from file
logger.info(f"Reading cell barcodes from {args.cell_file}")
cell_barcodes = set()
with open(args.cell_file, 'r') as f:
  for line in f:
    barcode = line.strip()
    if barcode:  # Skip empty lines
      cell_barcodes.add(barcode)

logger.info(f"Loaded {len(cell_barcodes)} cell barcodes")

# Extract reads matching the cell barcodes
logger.info(f"Reading from {args.input_bam}")
logger.info(f"Writing to {args.output_bam}")

reads_written = 0
reads_processed = 0

with pysam.AlignmentFile(args.input_bam, "rb") as inbam:
  header = inbam.header
  with pysam.AlignmentFile(args.output_bam, "wb", header=header) as outbam:
    for read in inbam.fetch(until_eof=True):
      reads_processed += 1
      
      if reads_processed % 1000000 == 0:
        logger.info(f"Processed {reads_processed:,} reads, written {reads_written:,} reads")
      
      # Check if read has the cell barcode tag
      if read.has_tag(args.cell_barcode_tag):
        barcode = read.get_tag(args.cell_barcode_tag)
        if barcode in cell_barcodes:
          outbam.write(read)
          reads_written += 1

logger.info(f"Finished processing {reads_processed:,} reads")
logger.info(f"Written {reads_written:,} reads to {args.output_bam}")

# Index the output BAM file if requested
if args.index:
  logger.info(f"Indexing {args.output_bam}")
  pysam.index(args.output_bam)
  logger.info("Indexing complete")