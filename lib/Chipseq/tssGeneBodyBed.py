import argparse
import sys
import logging
import os
import shutil
from gtfparse import read_gtf

def runCmd(cmd, logger):
  logger.info(cmd)
  os.system(cmd)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get bed file for TSS and gene body",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input gtf file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output prefix", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation.gtf"
  args.output = "/scratch/cqs_share/references/gencode/GRCh38.p13/gencode.v36.annotation"

logger = logging.getLogger('tssGeneBody')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("reading " + args.input + " ...")
with open(args.input, "rt") as fin:
  with open(args.output + ".tss.bed", "wt") as ftss:
    with open(args.output + ".genebody.bed", "wt") as fbody:
      for line in fin:
        if line.startswith("#"):
          continue

        parts = line.split('\t')
        if parts[2] != 'gene':
          continue
        
        infos = parts[8].split('; ')
        for info in infos:
          if info.startswith("gene_name"):
            gene_name = info[11:][:-1]
            tss_start = int(parts[3]) - 300
            tss_end = int(parts[3]) + 300
            ftss.write("%s\t%d\t%d\t%s\n" %(parts[0], tss_start, tss_end, gene_name))

            gb_start = tss_end + 1
            gb_end = int(parts[4]) + 3500
            fbody.write("%s\t%d\t%d\t%s\n" %(parts[0], gb_start, gb_end, gene_name))

logger.info("done.")