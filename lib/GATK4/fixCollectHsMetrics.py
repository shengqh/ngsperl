import argparse
import sys
import logging
import os
import random

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="fixCollectHsMetrics",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input CollectHsMetrics file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/PCA_scRNAseq/Exoseq/20220214_7538_CH/bwa_g4_refine_target_coverage/result/WD82458_NL_hs_metrics.txt"
  args.output = "/scratch/cqs/PCA_scRNAseq/Exoseq/20220214_7538_CH/bwa_g4_refine_target_coverage/result/WD82458_NL_hs_metrics.txt.fixed"

logger = logging.getLogger('fixCollectHsMetrics')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("reading " + args.input + " ...")

with open(args.output, "wt") as fout:
  with open(args.input, "rt") as fin:
    for line in fin:
      if "PCT_TARGET_BASES_100000X" not in line:
        fout.write(line)
      else:
        break
    
    parts = line.split('\t')
    PCT_TARGET_BASES_500X = parts.index('PCT_TARGET_BASES_500X')
    AT_DROPOUT = parts.index('AT_DROPOUT')
    del parts[PCT_TARGET_BASES_500X:AT_DROPOUT]
    fout.write("\t".join(parts))

    line = fin.readline()
    parts = line.split('\t')
    del parts[PCT_TARGET_BASES_500X:AT_DROPOUT]
    fout.write("\t".join(parts))

    for line in fin:
      fout.write(line)
  
logger.info("done.")
