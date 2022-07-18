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

def get_record_per_file(logger, input_files, trunk_number):
  total_line_count = 0
  for input_file in input_files:
    logger.info("Check read count of %s ..." % input_file)
    cur_line_count = rawgencount(input_file)
    logger.info(f"  {cur_line_count / 4} reads")
    total_line_count += cur_line_count
  total_record = total_line_count / 4
  result = math.ceil(total_record / trunk_number)
  logger.info("Total %d reads, each file should have almost %d reads." % (total_record, result))
  return (result)

def split_by_trunk(logger, inputFiles, isPairedEnd, outputFilePrefix, trunkNumber, startTrunk=1): 
  if isPairedEnd:
    read1files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 0]
    read2files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 1]
    input_array = [read1files, read2files]
  else:
    input_array = [inputFiles]

  recordPerFile = get_record_per_file(logger, input_array[0], trunkNumber)

  fileIndex = 0
  for read_files in input_array:
    fileIndex = fileIndex + 1
    trunk = startTrunk
    read_count = 0
    fout = None
    for inputFile in read_files:
      logger.info("Processing %s ..." % inputFile)
      
      with gzip.open(inputFile, "rt") as fin:
        while(True):
          line1 = fin.readline()
          if not line1:
            break
          line2 = fin.readline()
          line3 = fin.readline()
          line4 = fin.readline()

          if read_count == recordPerFile:
            fout.close()
            trunk += 1
            read_count = 0
          
          if read_count == 0:
            if isPairedEnd:
              fileName = f"{outputFilePrefix}.{str(trunk).zfill(2)}.{fileIndex}.fastq.gz"
            else:
              fileName = f"{outputFilePrefix}.{str(trunk).zfill(2)}.fastq.gz"
            logger.info("Writing reads to %s ..." % fileName)
            fout = gzip.open(fileName, "wt")

          read_count += 1
          fout.write(line1)
          fout.write(line2)
          fout.write(line3)
          fout.write(line4)
    fout.close()

  logger.info("done")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Scatter big FASTQ files to N small FASTQ files",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  #input can be multiple FASTQ files, such as F1_R1,F1_R2,F2_R1,F2_R2,F3_R1,F3_R2 in case some of them are more larger than others.
  #we don't know how many reads per sample
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', default="-", help="Output file prefix", required=NOT_DEBUG)
  parser.add_argument('--is_single_end', action='store', nargs='?', help="Is single end?")
  parser.add_argument('--trunk', action='store', nargs='?', type=int, default=50, help="Number of small files")
  
  args = parser.parse_args()
  
  if DEBUG:
    dfolder = "/scratch/vickers_lab/projects/20200708_smallRNA_KCV_3018_45_46_human_v5_byTiger/preprocessing/identical/result/"
    args.input = "%sCB10C_clipped_identical.fastq.gz,%sCB10C_clipped_identical.fastq.gz,%sCB11C_clipped_identical.fastq.gz,%sCB11C_clipped_identical.fastq.gz" % (dfolder, dfolder, dfolder, dfolder)
    args.outputPrefix = "/scratch/cqs/shengq2/temp/splitFastqTest"
    args.is_single_end = True
    args.trunk = 4
  
  logger = logging.getLogger('splitFastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  inputFiles = args.input.split(",")
  
  split_by_trunk(logger, inputFiles, not args.is_single_end, args.outputPrefix, args.trunk)
  
if __name__ == "__main__":
    main()
