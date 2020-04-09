import argparse
import sys
import logging
import os
import gzip
import math
from _ctypes import ArgumentError

def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)

def rawgencount(filename):
    with gzip.open(filename, 'rb') as f:
      f_gen = _make_gen(f.read)
      return sum( buf.count(b'\n') for buf in f_gen )

def remove_duplicate_paired_end(logger, read1file, read2file, outputFilePrefix):  
  logger.info("Processing %s ..." % read1file)
  totalLineCount = rawgencount(read1file)
  totalRecord = totalLineCount / 4
  logger.info("Total %d reads." % totalRecord)
  
  outfile1 = outputFilePrefix + ".1.fastq.gz"
  outfile2 = outputFilePrefix + ".2.fastq.gz"
  outfile1tmp = outfile1 + ".tmp.gz"
  outfile2tmp = outfile2 + ".tmp.gz"

  qnames = set()
  count = 0
  dupCount = 0
  with gzip.open(read1file, "rt") as fin1:
    with gzip.open(read2file, "rt") as fin2:
      with gzip.open(outfile1tmp, "wt") as fout1:
        with gzip.open(outfile2tmp, "wt") as fout2:
          while True:
            f1_1 = fin1.readline()
            if not f1_1:
              break

            f1_2 = fin1.readline()
            f1_3 = fin1.readline()
            f1_4 = fin1.readline()
            f2_1 = fin2.readline()
            f2_2 = fin2.readline()
            f2_3 = fin2.readline()
            f2_4 = fin2.readline()

            count += 1
            if count % 100000 == 0:
              logger.info("%d dup / %d processed / %d total ..." % (dupCount, count, totalRecord))

            qname = f1_1.rstrip().split(' ', 1)[0]
            if qname in qnames:
              logger.error("Repliacated query %s" % qname)
              dupCount += 1
              continue

            qnames.add(qname)

            fout1.write(f1_1)
            fout1.write(f1_2)
            fout1.write(f1_3)
            fout1.write(f1_4)
            fout2.write(f2_1)
            fout2.write(f2_2)
            fout2.write(f2_3)
            fout2.write(f2_4)

  os.rename(outfile1tmp, outfile1)        
  os.rename(outfile2tmp, outfile2) 
  logger.info("%d out of %d were replicated" % (dupCount, totalRecord))       
  logger.info("done")


def remove_duplicate_single_end(logger, read1file, outputFilePrefix):  
  logger.info("Processing %s ..." % read1file)
  totalLineCount = rawgencount(read1file)
  totalRecord = totalLineCount / 4
  logger.info("Total %d reads." % totalRecord)
  
  outfile1 = outputFilePrefix + ".fastq.gz"
  outfile1tmp = outfile1 + ".tmp.gz"

  qnames = set()
  count = 0
  dupCount = 0
  with gzip.open(read1file, "rt") as fin1:
    with gzip.open(outfile1tmp, "wt") as fout1:
      while True:
        f1_1 = fin1.readline()
        if not f1_1:
          break

        f1_2 = fin1.readline()
        f1_3 = fin1.readline()
        f1_4 = fin1.readline()

        count += 1
        if count % 100000 == 0:
          logger.info("%d / %d ..." % (count, totalRecord))

        qname = f1_1.rstrip().split(' ', 1)[0]
        if qname in qnames:
          logger.error("Repliacated query %s" % qname)
          dupCount += 1
          continue

        qnames.add(qname)

        fout1.write(f1_1)
        fout1.write(f1_2)
        fout1.write(f1_3)
        fout1.write(f1_4)

  os.rename(outfile1tmp, outfile1)        
  logger.info("%d out of %d were replicated" % (dupCount, totalRecord))       
  logger.info("done")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Split big fastq file to multiple small fastq files",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', default="-", help="Output file prefix", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/data/cqs/jennifer_pietenpol/20190905_3772/3772-LR-6-TCACTGTC-ATCTCCGG_S01_L001_R1_001.fastq.gz,/data/cqs/jennifer_pietenpol/20190905_3772/3772-LR-6-TCACTGTC-ATCTCCGG_S01_L001_R2_001.fastq.gz"
    args.outputPrefix = "/scratch/cqs/shengq2/temp/3772-LR-6"
  
  logger = logging.getLogger('dedupFastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  inputFiles = args.input.split(",")
  
  if len(inputFiles) == 1:
    remove_duplicate_single_end(logger, inputFiles[0], args.outputPrefix)
  elif len(inputFiles) == 2:
    remove_duplicate_paired_end(logger, inputFiles[0], inputFiles[1], args.outputPrefix)
  else:
    raise ArgumentError('inputFile should be only one file (single end) or two files (pair end): %s ' % args.input)
  
if __name__ == "__main__":
    main()
