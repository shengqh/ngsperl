import pysam
import argparse
import sys
import logging
from collections import OrderedDict

DEBUG = True
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Get read count in chromosomes",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM file', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output summary file", required=NOT_DEBUG)

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/jbrown_lab/shengq2/projects/20201208_chipseq_485_886_hg38/bowtie2_cleanbam/result/No_Treatment_886.noChrM.bam"
  args.output=args.input + ".chr_reads"

logger = logging.getLogger('bamReference')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

with pysam.Samfile(args.input, "rb") as sam:
  chr_map = OrderedDict()

  processed = 0
  for read in sam.fetch(until_eof=True):
    processed += 1
    if processed % 100000 == 0:
      logger.info(f"processed {processed}")

    if read.is_unmapped:
      continue

    if read.reference_name in chr_map:
      chr_map[read.reference_name] += 1
    else:
      chr_map[read.reference_name] = 1

with open(args.output, "wt") as fout:
  fout.write("Chromosome\tCount\n")
  for chr in chr_map.keys():
    fout.write("%s\t%d\n" % (chr, chr_map[chr]))
  
logger.info("done")