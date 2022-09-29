import argparse
import sys
import logging
import os
import os.path
import errno
import re
import gzip
from collections import OrderedDict

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get unmapped reads",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input fastq file', required=NotDEBUG)
parser.add_argument('-m', '--mapped', action='store', nargs='?', help="Input mapped read count list file", required=NotDEBUG)
parser.add_argument('-c', '--count', action='store', nargs='?', help="Input count file", required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output unmapped FASTQ file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/vickers_lab/projects/20220830_8643_CM_smRNA_human_bakeoff/intermediate_data/bowtie1_genome_unmapped_reads/result/Pool1_Exomere_1_clipped_identical.unmapped.fastq.gz"
  args.mapped = "/scratch/vickers_lab/projects/20220830_8643_CM_smRNA_human_bakeoff/final_unmapped/final_unmapped_reads_python/result/CM_8643_bakeoff__fileList2.list"
  args.output = "/scratch/vickers_lab/projects/20220830_8643_CM_smRNA_human_bakeoff/final_unmapped/final_unmapped_reads_python/result/Pool1_Exomere_1.unmapped.fastq.gz"

logger = logging.getLogger('unmappedReads')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

check_file_exists(args.input)
check_file_exists(args.mapped)
check_file_exists(args.count)

mapped = set()
with open(args.mapped, "rt") as fin:
  for line in fin:
    read_file = line.split('\t')[0]
    logger.info(f"reading {read_file}")
    with open(read_file, "rt") as fr:
      fr.readline()
      for seq in fr:
        read = seq.split('\t', 1)[0]
        if read not in mapped:
          mapped.add(read)

logger.info("writing unmapped reads to " + args.output + " ...")
with gzip.open(args.output, "wt") as fout:
  with gzip.open(args.input, "rt") as fin:
    icount = 0
    while True:
      header = fin.readline()
      if '' == header:
        break

      if not header.startswith("@"):
        continue

      icount += 1
      if icount % 100000 == 0:
        logger.info(icount)

      sequence = fin.readline().rstrip()
      line3 = fin.readline()
      line4 = fin.readline()

      if sequence in mapped:
        continue

      fout.write(header)
      fout.write(f"{sequence}\n")
      fout.write(line3)
      fout.write(line4)

dupcount = re.sub(".gz$", ".dupcount", args.output)
logger.info("writing unmapped reads to " + dupcount + " ...")
with open(dupcount, "wt") as fout:
  with open(args.count, "rt") as fin:
    fout.write(fin.readline())
    icount = 0
    for line in fin:
      icount += 1
      if icount % 100000 == 0:
        logger.info(icount)

      parts = line.rstrip().split('\t')
      if parts[2] in mapped:
        continue
      fout.write(line)

logger.info("done")