import argparse
import sys
import logging
import os

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Remove chromosome with only 1 segment",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input interval file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output interval file", required=NotDEBUG)
parser.add_argument('-c', '--contig_ploidy_priors_file', action='store', nargs='?', help="Contig ploidy priors file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_03_FilterIntervals/result/human_exomeseq.filtered.all.interval_list"
  args.output = "/scratch/weissvl/shengq2/20210616_human_exomeseq/bwa_g4_refine_gatk4_CNV_Germline_03_FilterIntervals/result/human_exomeseq.filtered.interval_list"
  args.contig_ploidy_priors_file = "/data/cqs/references/broad/hg38/contig_ploidy_priors_homo_sapiens.chr.tsv"

logger = logging.getLogger('filterIntervals')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

chromSet = set()
with open(args.contig_ploidy_priors_file, "rt") as fin:
  fin.readline()
  for line in fin:
    parts = line.split('\t')
    chromSet.add(parts[0])

with open(args.output, "wt") as fout:
  with open(args.input) as fh:
    for line in fh:
      if line.startswith("@"):
        if line.startswith("@SQ"):
          parts = line.split('\t')
          chrom = parts[1].split(':')[1]
          if chrom in chromSet:
            fout.write(line)
        else:
          fout.write(line)
        continue
      
      parts = line.split('\t', 1)
      if parts[0] in chromSet:
        fout.write(line)

logger.info("done.")