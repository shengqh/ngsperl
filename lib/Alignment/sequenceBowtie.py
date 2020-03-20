import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from BowtieIndex import BowtieIndexItem
from _ctypes import ArgumentError

def readBowtieList(bowtieIndexList):
  result = []
  with open(bowtieIndexList) as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      if len(parts) >= 2:
        result.append(BowtieIndexItem(parts[0], parts[1]))
  return(result)

def alignment(logger, inputFile, outputFilePrefix, bowtieIndexList, threadNumber):
  bowtieIndexList = readBowtieList(bowtieIndexList)

  categories = sorted(set([bi.Category for bi in bowtieIndexList])) 
  logfile = outputFilePrefix + ".log"
  with open(outputFilePrefix + ".txt", "wt") as fout:
    with open(logfile, "wt") as flog:
      for category in categories:
        logger.info("Searching category %s ..." % category)
        flog.write(">" + category + "\n")
        flog.flush()
        bowtieIndecies = [bi.Index for bi in bowtieIndexList if bi.Category == category]

        bowtieCount = 0

        inputFastq = inputFile
        for bowtieIndex in bowtieIndecies:
          bowtieCount = bowtieCount + 1
          logger.info("  Searching to %s ..." % bowtieIndex)
          outfile = "%s.%d.out" % (outputFilePrefix, bowtieCount)
          unmapped = "%s.%d.unmapped.fastq" % (outputFilePrefix, bowtieCount)

          subprocess.call(['bowtie', '-k', '1', '-v', '0', '-p', str(threadNumber), '--no-unal', '--un', unmapped, bowtieIndex, inputFastq, outfile], 
            stderr=flog)

          with open(outfile, "rt") as fin:
            for line in fin:
              parts = line.rstrip().split('\t')
              fout.write("%s\t%s\t%s\t%s\t%s\n" % (parts[0].split(' ')[0], parts[1], parts[2], parts[3], category))

          os.remove(outfile)
          
          if bowtieCount > 1:
            os.remove(inputFastq)
          inputFastq = unmapped
        
        if inputFastq != inputFile:
          os.remove(inputFastq)
  
  logger.info("done")

def main():
  DEBUG = True
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Align fastq file to a set of bowtie index and combine the result",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-d', '--databaseFileList', action='store', nargs='?', help='Input file list of bowtie indecies', required=NOT_DEBUG)
  parser.add_argument('-t', '--threadNumber', action='store', nargs='?', type=int, default=8, help="Thread number")
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', default="-", help="Output file prefix", required=NOT_DEBUG)
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_clipped_identical.unmapped.fastq.gz"
    args.outputPrefix = "/scratch/cqs/shengq2/temp/VLDL_WZ_bacteria"
    args.databaseFileList = "/scratch/cqs/kasey_vickers_projects/testdata/bacteria.list"
  
  logger = logging.getLogger('sequenceBowtie')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  alignment(logger, args.input, args.outputPrefix, args.databaseFileList, args.threadNumber)
  
if __name__ == "__main__":
    main()
