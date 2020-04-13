import argparse
import sys
import logging
import os
import math

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Generate scatter interval files",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input interval file', required=NotDEBUG)
parser.add_argument('-n', '--scatter_number', type=int, default=100, help="Input scatter number", required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file prefix", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/cqs/shengq2/jennifer/20200407_lindsay_exomeseq_3772_hg38/bwa_refine_nosoftclip_gatk4_CNV_Germline_03_FilterIntervals/result/lindsay_exomeseq_3772.filtered.interval_list"
  args.output = "/scratch/cqs/shengq2/jennifer/20200407_lindsay_exomeseq_3772_hg38/bwa_refine_nosoftclip_gatk4_CNV_Germline_03_FilterIntervals/result/scatter/lindsay_exomeseq_3772"

logger = logging.getLogger('scatterInterval')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

headers = []
lines = []
with open(args.input, "rt") as fin:
  for line in fin:
    if line.startswith("@"):
      headers.append(line)
    else:
      lines.append(line)

maxLength = len(str(args.scatter_number))
numberEachFile = math.ceil(len(lines) / args.scatter_number)
for idx in range(1, (args.scatter_number+1)):
  startIndex = (idx-1) * numberEachFile
  endIndex = idx * numberEachFile
  endIndex = min(endIndex, len(lines))

  #idxstr = str(idx).zfill(maxLength)
  idxstr = str(idx)
  subfile = args.output + "." + idxstr + ".interval_list"
  #print("%s:%d-%d" %(idxstr, startIndex, endIndex))
  with open(subfile, "wt") as fout:
    for header in headers:
      fout.write(header)
    for rindex in range(startIndex, endIndex):
      fout.write(lines[rindex])

logger.info("done")