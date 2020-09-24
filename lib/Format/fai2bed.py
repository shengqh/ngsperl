import pysam
import argparse
import sys
import logging
import os
from asyncore import read

parser = argparse.ArgumentParser(description="Convert fai to bed",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fai file', required=True)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output bed file", required=True)

args = parser.parse_args()

logger = logging.getLogger('fai2bed')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("writing to %s" % args.output)

with open(args.input, "rt") as fin:
  with open(args.output, "wt") as fout:
    for line in fin:
      parts = line.strip().split('\t')
      chrom = parts[0]
      chromLength = parts[1]
      fout.write(f"{chrom}\t0\t{chromLength}\n")

logger.info("done")
