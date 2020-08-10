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

def check_fastq_duplicate(logger, inputFile, outputFilePrefix): 
  inputFiles = inputFile.split(',')
  inputFile = inputFiles[0]

  logger.info("Processing %s ..." % inputFile)
  ignore = False
  querySet = {}
  dupQuery = {}
  index = 0
  dupCount = 0
  with gzip.open(inputFile, "rt") as fin:
    while(True):
      line1 = fin.readline()
      if not line1:
        break
      fin.readline()
      fin.readline()
      fin.readline()
      index += 1

      if index % 1000000 == 0:
        logger.info(index)

      parts = line1.split(' ', 1)
      id = parts[0].rstrip()
      if id in querySet:
        dupCount += 0

        if ignore:
          continue

        if id in dupQuery:
          dupQuery[id].append(index)
        else:
          dupQuery[id] = [querySet[id], index]

        if len(dupQuery) == 10:
          logger.info("Detect 10+ duplicates")
          ignore = True
      else:
        querySet[id] = index
        
  with open(outputFilePrefix + ".txt", 'wt') as fout:
    fout.write("Result\tTotal\tDuplicated\n")
    result = "FAIL" if len(dupQuery) > 0 else "PASS"
    fout.write("%s\t%d\t%d\n" % (result, index, dupCount))

  if len(dupQuery) > 0:
    first10file = outputFilePrefix + ".first10.txt"
    with open(first10file, "wt") as fout:
      fout.write("Query\tPosition\n")
      dqs = sorted(dupQuery.items(), key=lambda x:x[1][0])
      logger.info("There are duplicated reads. First 10 were saved in %s" % first10file)
      for dq in dqs:
        id = dq[0]
        queries = dq[1]
        fout.write("%s\t%s\n" % (id, ",".join(str(q) for q in queries)))
      
  logger.info("done")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Check duplicates in fastq file.",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq file', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', help="Output summary file prefix", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/data/cqs/ramirema/ciombor_kristen_data/Plate3/2585-KL-215-GCAATATT-GACTGAGT_S51_R1_001.fastq.gz"
    args.output = "/scratch/cqs/shengq2/temp/duplicate"
  
  logger = logging.getLogger('check_fastq_duplicate')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  check_fastq_duplicate(logger, args.input, args.outputPrefix)
  
if __name__ == "__main__":
    main()
