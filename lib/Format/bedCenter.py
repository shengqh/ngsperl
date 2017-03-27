import pysam
import argparse
import sys
import logging
import os
from asyncore import read

parser = argparse.ArgumentParser(description="Center bed with assigned window.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG = True
NOT_DEBUG = not DEBUG

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bed file (use "-" as stdin)', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output bed file", required=NOT_DEBUG)
parser.add_argument('-w', '--window', action='store', nargs='?', type=int, default=5000, help="Window at each side of peak center")

args = parser.parse_args()

if DEBUG:
  args.input="/scratch/cqs/shengq1/brown/20170317_chipseq_gse53998_hg19/macs1callpeak/result/EC_BRD4_CON/EC_BRD4_CON_peaks.name.bed"
  args.window=5000
  args.output="/scratch/cqs/shengq1/brown/20170317_chipseq_gse53998_hg19/macs1callpeak/result/EC_BRD4_CON_peaks.center.bed"

logger = logging.getLogger('bedCenter')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

tmpfile = args.output + ".tmp"
with open(tmpfile, 'w') as out:
  with open(args.input, "r") as ins:
    for line in ins:
      parts = line.strip().split('\t')
      center = (int(parts[1]) + int(parts[2])) / 2
      left=center - args.window
      if left >= 0:
        out.write("%s\t%d\t%d\t%s\n" % (parts[0], left, center + args.window, "\t".join(parts[3:])))
  
if os.path.isfile(args.output):
  os.remove(args.output)
os.rename(tmpfile, args.output)
        
  
