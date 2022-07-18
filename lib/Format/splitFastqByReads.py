import argparse
import sys
import logging
import os
import gzip
import math
from _ctypes import ArgumentError

def split(logger, inputFiles, is_single_end, outputFilePrefix, recordPerFile): 
  if is_single_end:
    input_array = [inputFiles]
  else:
    read1files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 0]
    read2files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 1]
    input_array = [read1files, read2files]

  fileIndex = 0
  for read_files in input_array:
    fileIndex = fileIndex + 1
    trunk = 1
    read_count = 0
    fout = None
    for inputFile in read_files:
      logger.info("Processing %s ..." % inputFile)
      
      with gzip.open(inputFile, "rt") as fin:
        while(True):
          line1 = fin.readline()
          if not line1:
            break

          if line1[0] != '@':
            raise Exception("Wrong id: %s" % line1)
          
          line2 = fin.readline()
          line3 = fin.readline()
          line4 = fin.readline()

          if read_count == recordPerFile:
            fout.close()
            trunk += 1
            read_count = 0
          
          if read_count == 0:
            if isPairedEnd:
              fileName = "%s.%d.%d.fastq.gz" % (outputFilePrefix, trunk, fileIndex)
            else:
              fileName = "%s.%d.fastq.gz" % (outputFilePrefix, trunk)
            logger.info("Writing reads to %s ..." % fileName)
            fout = gzip.open(fileName, "wt")

          read_count += 1
          fout.write(line1)
          fout.write(line2)
          fout.write(line3)
          fout.write(line4)
    if fount != None:
      fout.close()

  logger.info("done")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Split big fastq files to fixed number of reads per small FASTQ file",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  #We don't know how many files would be generated
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', default="-", help="Output file prefix", required=NOT_DEBUG)
  parser.add_argument('--is_single_end', action='store', nargs='?', help="Is single end?")
  parser.add_argument('--record_per_file', action='store', nargs='?', type=int, default=5000000, help="Record per file")
  
  args = parser.parse_args()
  
  if DEBUG:
    dfolder = "/scratch/vickers_lab/projects/20200708_smallRNA_KCV_3018_45_46_human_v5_byTiger/preprocessing/identical/result/"
    args.input = "%sCB10C_clipped_identical.fastq.gz,%sCB10C_clipped_identical.fastq.gz,%sCB11C_clipped_identical.fastq.gz,%sCB11C_clipped_identical.fastq.gz" % (dfolder, dfolder, dfolder, dfolder)
    args.outputPrefix = "/scratch/cqs/shengq2/temp/splitFastqTest"
    args.is_single_end = True
    args.record_per_file = 5000000
  
  logger = logging.getLogger('splitFastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  inputFiles = args.input.split(",")
  
  split(logger, inputFiles, args.is_single_end, args.outputPrefix, args.record_per_file)
  
if __name__ == "__main__":
    main()
