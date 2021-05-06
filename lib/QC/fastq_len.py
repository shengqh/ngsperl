import argparse
import sys
import logging
import os
import gzip
from asyncore import read

parser = argparse.ArgumentParser(description="Get length distribution in FASTQ file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

DEBUG=False
NOT_DEBUG=not DEBUG

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input FASTQ files (use "," as seperation)', required=NOT_DEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file name")
args = parser.parse_args()

if DEBUG:
  args.input = "/scratch/stein_lab/20210429_qiong_rnaseq_6130_gut_hg38/cutadapt/result/S324_DC_RA_RNA_clipped.1.fastq.gz,/scratch/stein_lab/20210429_qiong_rnaseq_6130_gut_hg38/cutadapt/result/S324_DC_RA_RNA_clipped.2.fastq.gz"
  args.output = "/scratch/stein_lab/20210429_qiong_rnaseq_6130_gut_hg38/cutadapt/test.len"

logger = logging.getLogger('fastq_len')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

fastq_files = args.input.split(',')

result={}
for fq in fastq_files:
  logger.info("Processing %s ..." % fq)
  reads_count = 0
  fin = gzip.open(fq, "rt") if fq.endswith(".gz") else open(fq, "rt")
  with fin:
    while True:
      query = fin.readline()
      if not query:
        break
      if not query.startswith("@"):
        raise Exception("Wrong format: " + query)

      reads_count += 1
      if reads_count % 1000000 == 0:
        logger.info("%d reads processed." % reads_count)
        #break

      sequence = fin.readline().rstrip()
      fin.readline()
      fin.readline()

      seqlen = len(sequence)
      result[seqlen] = result.setdefault(seqlen, 0) + 1

with open(args.output, "wt") as fout:
  fout.write("Len\tCount\n")
  maxlen = max(result.keys())
  for len in range(1, maxlen):
    fout.write("%d\t%d\n" % (len, result.setdefault(len, 0)))

logger.info("done")