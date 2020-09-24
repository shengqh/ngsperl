import sys
import gzip
import os
import logging
import argparse
import subprocess

DEBUG = False
NOT_DEBUG= not DEBUG

parser = argparse.ArgumentParser(description="Get first N reads from FASTQ file.",
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input FASEQ file (use "," to join pair end files)', required=NOT_DEBUG)
parser.add_argument('-p', '--pair_end', action='store_true', help='Is pair end data?')
parser.add_argument('-n', '--number', action='store', type=int, nargs='?', default=10000, help='Input number of reads')
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = ""
  args.output = ""
  
print(args)

sourceFile = args.input
isPairEnd = args.pair_end
number = args.number
outputFile = args.output

logger = logging.getLogger('getFirstNReads')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def extract(sourceFile, extractCount, targetFile):
  linecount = extractCount * 4
  with gzip.open(sourceFile, "rt") as fin:
    with gzip.open(targetFile, "wt") as fout:
      for count in range(0, linecount):
        fout.write(fin.readline())

if isPairEnd:
  f1 = sourceFile.split(',')[0]
  f2 = sourceFile.split(',')[1]
  read1 = outputFile.replace("2.fastq.gz", "1.fastq.gz")
  extract(f1, number, read1)
  extract(f2, number, outputFile)
else:
  extract(sourceFile, number, outputFile)

logger.info("done.")
