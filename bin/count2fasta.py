#!/usr/bin/env python3
 
import os
import os.path
import argparse
import sys
import logging
import re
  
DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Convert count table to fasta.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input count table file', required=NotDEBUG)
parser.add_argument('-m', '--minimum_count', type=int, action='store', nargs='?', default=20, help="Minimum count required")
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output fasta file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/workspace/shengq2/kasey_vickers/2019_projects/20190905_3791_RA2_Run2_mouse_V4T/host_genome/bowtie1_genome_host_reads_table/result/RA2_3791_Run2.count"
  args.output = "/workspace/shengq2/kasey_vickers/2019_projects/20190905_3791_RA2_Run2_mouse_V4T/host_genome/bowtie1_genome_host_reads_table/result/RA2_3791_Run2.fasta"

logger = logging.getLogger('count2fasta')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

index = 0
discard = 0

saved=[]
with open(args.input, "rt") as fin:
  line = fin.readline()
  for line in fin:
    parts = line.rstrip().split('\t')
    total_count = sum(float(p) for p in parts[1:])
    if total_count < args.minimum_count:
      discard = discard + 1
      continue

    saved.append([parts[0], total_count])

saved = sorted(saved, key=lambda read: read[1], reverse=True)
with open(args.output, "wt") as fout:
  for idx, read in enumerate(saved):
    fout.write(">%d_%d\n%s\n" % (idx+1, read[1], read[0]))

logger.info("Saved=%d, discarded=%d" % (len(saved), discard))
logger.info("done.")
