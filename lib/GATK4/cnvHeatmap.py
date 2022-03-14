import argparse
import sys
import logging
import os
import random

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="CNV heatmap",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input combined CNV files', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/jennifer/20190906_lindsay_exomeseq_3772_hg38/GATK4_CNV_Germline_07_CombineGCNV/result/lindsay_exomeseq_3772.txt"
  args.output = "/scratch/cqs/shengq2/jennifer/20190906_lindsay_exomeseq_3772_hg38/GATK4_CNV_Germline_07_CombineGCNV/result/lindsay_exomeseq_3772.heatmap.txt"

logger = logging.getLogger('cnvHeatmap')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

random.seed(20190930)

logger.info("reading " + args.input + " ...")

with open(args.input, "rt") as fin:
  headers = fin.readline().split('\t')
  samples = headers[2:]
  with open(args.output, "wt") as fout:
    fout.write("chr\tstart\tend\tmedian.bp\t" + "\t".join(samples))
    for line in fin:
      parts = line.split('\t')

      list_set = set(part.split(',')[0] for part in parts[5:])
      if len(list_set) == 1:
        continue

      chrom = parts[0]
      start = int(parts[1])
      end = int(parts[2])
      median = (start + end) / 2
      fout.write("%s\t%d\t%d\t%d" % (chr, start, end, median))
      for part in parts[5:]:
        part = part.rstrip()
        uf = random.uniform(-0.01, 0.01)
        if part == "":
          fout.write("\t%.4f" % uf)
          continue

        cnvparts = part.split(',')

        if cnvparts[2] == '0':
          fout.write("\t%.4f" % (uf - 1))
          continue

        if cnvparts[2] == '1':
          fout.write("\t%.4f" % (uf - 0.5))
          continue

        fout.write("\t%.1f" % (uf + float(cnvparts[2]) / 2.0))
      fout.write('\n')
  
logger.info("done.")
