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

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/jennifer/20200407_lindsay_exomeseq_3772_hg38/bwa_refine_nosoftclip_gatk4_CNV_Germline_03_FilterIntervals/result/lindsay_exomeseq_3772.filtered.interval_list.old"
  args.output = "/scratch/cqs/shengq2/jennifer/20200407_lindsay_exomeseq_3772_hg38/bwa_refine_nosoftclip_gatk4_CNV_Germline_03_FilterIntervals/result/lindsay_exomeseq_3772.filtered.interval_list.old.txt"

logger = logging.getLogger('filterIntervals')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

lastchrom = ""
chroms = []
with open(args.output, "wt") as fout:
  with open(args.input) as fh:
    for line in fh:
      if line.startswith("@"):
        fout.write(line)
        continue
      
      parts = line.split('\t', 1)
      if parts[0] != lastchrom:
        if len(chroms) > 1:
          for bed in chroms:
            fout.write(bed)
        chroms = [line]
        lastchrom = parts[0]
      else:
        chroms.append(line)
  
  if len(chroms) > 1:
    for bed in chroms:
      fout.write(bed)

logger.info("done.")