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

def split_by_trunk(logger, inputFiles, isPairedEnd, outputFilePrefix, trunkNumber, startTrunk=1, fill_length=3, compresslevel=1, total_read=0): 
  if isPairedEnd:
    read1files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 0]
    read2files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 1]
    input_array = [read1files, read2files]
  else:
    input_array = [inputFiles]

  if total_read > 0:
    recordPerFile = math.ceil(total_read / trunkNumber)
  else:
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
              fileName = f"{outputFilePrefix}.{str(trunk).zfill(fill_length)}.{fileIndex}.fastq.gz"
            else:
              fileName = f"{outputFilePrefix}.{str(trunk).zfill(fill_length)}.fastq.gz"
            logger.info("Writing reads to %s ..." % fileName)
            fout = gzip.open(fileName, "wt",  compresslevel=compresslevel)

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
  parser.add_argument('--total_reads', action='store', nargs='?', type=int, default=0, help="Total reads in input")
  parser.add_argument('--start_trunk', action='store', nargs='?', type=int, default=1, help="Trunk starts from (default 1)")
  parser.add_argument('--fill_length', action='store', nargs='?', type=int, default=3, help="Trunk name length (fill with zero)")
  parser.add_argument('--compresslevel', action='store', nargs='?', type=int, default=1, help="Compress level, 1: fastest, 9: slowest")

  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/data/jbrown_lab/2023/20230602_10026_wgs_Human_Cholangiocyte_Cell/10026-DB-0001_S1_L005_R1_001.fastq.gz,/data/jbrown_lab/2023/20230602_10026_wgs_Human_Cholangiocyte_Cell/10026-DB-0001_S1_L005_R2_001.fastq.gz"
    args.outputPrefix = "/nobackup/brown_lab/projects/20230602_wgs_10026_hg38/bwa_00_splitFastq/result/DB_0001"
    args.is_single_end = False
    args.total_reads = 1168866345
    args.trunk = 20
  
  logger = logging.getLogger('splitFastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  inputFiles = args.input.split(",")
  
  split_by_trunk(
    logger = logger, 
    inputFiles = inputFiles,
    isPairedEnd= not args.is_single_end,
    outputFilePrefix= args.outputPrefix, 
    trunkNumber= args.trunk,
    startTrunk=args.start_trunk,
    fill_length= args.fill_length,
    compresslevel= args.compresslevel,
    total_read= args.total_reads)
  
if __name__ == "__main__":
    main()
