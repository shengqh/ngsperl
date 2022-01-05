import sys
import gzip
import os
import logging
import argparse
import subprocess

DEBUG = False
NOT_DEBUG= not DEBUG

if DEBUG:
  sourceFile="/scratch/jbrown_lab/shengq2/projects/20200402_chipseq_4615_human_cqs/test_extract/result/Chipseq_4615_liver__fileList1.list"
  isPairEnd=True
  outputFile="/scratch/jbrown_lab/shengq2/projects/20200402_chipseq_4615_human_cqs/test_extract/result/Chipseq_4615_liver.2.fastq.gz"
else:
  parser = argparse.ArgumentParser(description="Generate smallRNA count from count xml.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq list file', required=NOT_DEBUG)
  parser.add_argument('-p', '--pair_end', action='store_true', help='Is paired end?')
  parser.add_argument('-o', '--output', action='store', nargs='?', help="Output Fastq file", required=NOT_DEBUG)

  args = parser.parse_args()
  
  print(args)
  
  sourceFile = args.input
  isPairEnd = args.pair_end
  outputFile = args.output

logger = logging.getLogger('extractFirstNFastq')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

def getFirstFileList(fileName, isPairEnd):
  result = None
  with open(fileName) as fh:
    for line in fh:
      read1 = line.strip().split('\t', 1)[0]
      if isPairEnd:
        read2 = line.strip().split('\t', 1)[0]
        result = [read1, read2]
      else:
        result = [read1]
      break
  return(result)

def extract(sourceFile, targetFile):
  extractCount = 10000
  with gzip.open(sourceFile, "rt") as fin:
    with gzip.open(targetFile, "wt") as fout:
      for count in range(0, extractCount * 4):
        fout.write(fin.readline())

fastqFiles = getFirstFileList(sourceFile, isPairEnd)
if isPairEnd:
  read1 = outputFile
  read2 = read1.replace("1.fastq.gz", "2.fastq.gz")
  extract(fastqFiles[0], read1)
  extract(fastqFiles[1], read2)
else:
  extract(fastqFiles[0], outputFile)

logger.info("done.")
