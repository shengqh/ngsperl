import subprocess
import os.path
import re
import argparse
import sys
from pybedtools import BedTool

DEBUG = True

parser = argparse.ArgumentParser(description="find overlap gene.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

if not DEBUG:
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input locus file (bed format)', required=True)
  parser.add_argument('-g', '--gene_file', action='store', nargs='?', help='Gene locus file (bed format)', required=True)
  parser.add_argument('-o', '--output', action='store', nargs='?', help='Output overlap file', required=True)

  args = parser.parse_args()
  input_file=args.input
  gene_file = args.gene_file
  output_file=args.output
else:
  input_file= "/scratch/cqs/shengq1/vickers/20170720_AGO_human_CLIP/macs2/result/GSM1020022/GSM1020022_peaks.narrowPeak.bed"
  gene_file = "/scratch/cqs/shengq1/references/smallrna/v3/hg19_miRBase21_GtRNAdb2_gencode19_ncbi.sorted.bed"
  output_file="/scratch/cqs/shengq1/vickers/20170720_AGO_human_CLIP/macs2/result/GSM1020022/GSM1020022_peaks.narrowPeak.overlap.tsv"

closet = [nearest for nearest in BedTool(input_file).closest(gene_file, d=True)]

with open(output_file, 'w') as w:
  for nearest in closet:
    overlap = nearest.fields[12]
    if overlap == u'0':
      w.write(str(nearest))
