import argparse
import sys
import logging
import os
from Bio import SeqIO
import subprocess

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Find duplicated sample",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input plink smiss file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output plink duplicated sample list file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/macrae_linton/20191008_linton_megachip_3778_human/T01_plinkqc/result/linton_rmdup_snp.smiss"
  args.output = "/scratch/cqs/shengq2/macrae_linton/20191008_linton_megachip_3778_human/T01_plinkqc/result/linton_rmdup_snp.smiss.txt"

logger = logging.getLogger('findDuplicatedSample')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

idmap = {}
with open(args.input, "rt") as fin:
  header = fin.readline()
  for line in fin:
    parts = line.split('\t')
    if not parts[1] in idmap:
      idmap[parts[1]] = [[float(parts[4]), parts[0], parts[1]]]
    else:
      idmap[parts[1]].append([float(parts[4]), parts[0], parts[1]])

with open(args.output, "wt") as fout:
  for iid in idmap.keys():
    values = idmap[iid]
    if len(values) == 1:
      continue

    minIndex = 0
    minValue = values[0][0]
    for idx in range(len(values)):
      vi = values[idx]
      if vi[0] < minValue:
        minIndex = idx
        minValue = vi[0]
    
    for idx in range(len(values)):
      if idx != minIndex:
        vi = values[idx]
        fout.write("%s\t%s\n" % (vi[1], vi[2]))

logger.info("Done.")
