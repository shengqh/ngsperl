import argparse
import sys
import logging
import os
import os.path
import errno
import re
import gzip
from collections import OrderedDict
from CountXmlUtils import readCountXmlQueryNames
from DupCountUtils import readDupCountMap

def check_file_exists(file):
  if not os.path.exists(file):
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Get unmapped reads",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input bam file', required=NotDEBUG)
parser.add_argument('-s', '--smRNA', action='store', nargs='?', help="Input smallRNA count xml file", required=NotDEBUG)
parser.add_argument('-p', '--pmHost', action='store', nargs='?', help="Input perfect mapped name file", required=NotDEBUG)
parser.add_argument('-c', '--count', action='store', nargs='?', help="Input dupcount file", required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output unmapped FASTQ file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/vickers_lab/projects/20201201_5188_ES_smRNA_hg38_v5_byMars/preprocessing/identical/result/HDL_96_clipped_identical.fastq.gz"
  args.smRNA = "/scratch/vickers_lab/projects/20201201_5188_ES_smRNA_hg38_v5_byMars/intermediate_data/bowtie1_genome_1mm_NTA_smallRNA_count/result/HDL_96/HDL_96.count.mapped.xml"
  args.pmHost = "/scratch/vickers_lab/projects/20201201_5188_ES_smRNA_hg38_v5_byMars/intermediate_data/bowtie1_genome_1mm_NTA_pmnames/result/HDL_96.pmnames"
  args.count = "/scratch/vickers_lab/projects/20201201_5188_ES_smRNA_hg38_v5_byMars/preprocessing/identical/result/HDL_96_clipped_identical.fastq.dupcount"
  args.output = "/scratch/vickers_lab/projects/20201203_tellseq_hg38/HDL_96.unmapped.fastq.gz"

logger = logging.getLogger('unmappedReads')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

check_file_exists(args.input)
check_file_exists(args.smRNA)
check_file_exists(args.pmHost)
check_file_exists(args.count)

logger.info("reading " + args.count + " ...")
countDict = readDupCountMap(args.count)

logger.info("reading " + args.smRNA + " ...")
smRNAReads = readCountXmlQueryNames(args.smRNA)
smRNAqueries = set([sr.split(":CLIP_")[0] for sr in smRNAReads])

logger.info("reading " + args.pmHost + " ...")
pmHosts = set()
with open(args.pmHost, "rt") as fin:
  fin.readline()
  for line in fin:
    qname = line.split(':CLIP_')[0]
    if qname not in smRNAqueries:
      pmHosts.add(qname)

mappedReads = smRNAqueries.union(pmHosts)

unmappedCount = 0
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

      sequence = fin.readline()
      line3 = fin.readline()
      line4 = fin.readline()

      id = header.rstrip()[1:].split(' ')[0].split('\t')[0]

      if id in mappedReads:
        continue

      unmappedCount += countDict[id]

      fout.write(header)
      fout.write(sequence)
      fout.write(line3)
      fout.write(line4)

with open(args.output + '.info', 'wt') as fout:
  fout.write("Category\tCount\n")
  fout.write("smallRNA\t%d\n" % sum([countDict[r] for r in smRNAqueries]))
  fout.write("pmHost\t%d\n" % sum([countDict[r] for r in pmHosts]))
  fout.write("unmapped\t%d\n" % unmappedCount)

logger.info("done")