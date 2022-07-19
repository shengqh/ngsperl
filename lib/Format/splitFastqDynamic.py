import argparse
import logging
import os
import math
import enum
import subprocess
from splitFastq import split_by_trunk

# Enum for size units
class SIZE_UNIT(enum.Enum):
   BYTES = 1
   KB = 2
   MB = 3
   GB = 4

def convert_unit(size_in_bytes, unit):
   """ Convert the size from bytes to other units like KB, MB or GB"""
   if unit == SIZE_UNIT.KB:
       return size_in_bytes/1024
   elif unit == SIZE_UNIT.MB:
       return size_in_bytes/(1024*1024)
   elif unit == SIZE_UNIT.GB:
       return size_in_bytes/(1024*1024*1024)
   else:
       return size_in_bytes

def get_file_size(file_name, size_type = SIZE_UNIT.GB ):
   """ Get file in size in given unit like KB, MB or GB"""
   size = os.path.getsize(file_name)
   return convert_unit(size, size_type)

def do_symlink(source, target):
  if os.path.exists(target):
    os.remove(target)
  
  os.symlink(source, target)

def do_fastqsplitter(logger, input_fastq, start_trunk, expect_trunk, output_file, fill_length, compresslevel=1):
  options = ['fastqsplitter', '-i', input_fastq, "-c", compresslevel]
  for trunk in range(start_trunk, start_trunk + expect_trunk):
    trunk_str = str(trunk).zfill(fill_length)
    trunk_file = output_file.replace("__TRUNK__", trunk_str)
    options = options + ["-o", trunk_file]
  logger.info(str(options))
  subprocess.call(options)

def split_dynamic_paired_end(logger, inputFiles, outputFilePrefix, min_file_size_gb, trunk_file_size_gb, fill_length, compresslevel=1, call_fastqsplitter=False): 
  input_array = [[inputFiles[idx], inputFiles[idx+1]] for idx in range(0, len(inputFiles), 2)]

  trunk = 1
  for read_files in input_array:
    read1 = read_files[0]
    read2 = read_files[1]
    read1_gb = get_file_size(read1)

    if read1_gb < min_file_size_gb:
      file1 = f"{outputFilePrefix}.{str(trunk).zfill(fill_length)}.1.fastq.gz"
      file2 = f"{outputFilePrefix}.{str(trunk).zfill(fill_length)}.2.fastq.gz"
      logger.info(f"softlink {read1} => {file1}")
      logger.info(f"softlink {read2} => {file2}")
      do_symlink(read1, file1)
      do_symlink(read2, file2)
      trunk += 1
      continue

    expect_trunk = math.ceil(read1_gb / trunk_file_size_gb)
    if not call_fastqsplitter:
      split_by_trunk(logger, read_files, True, outputFilePrefix, expect_trunk, trunk, fill_length, compresslevel)
    else:
      do_fastqsplitter(logger, read1, trunk, expect_trunk, f"{outputFilePrefix}.__TRUNK__.1.fastq.gz", fill_length, compresslevel)
      do_fastqsplitter(logger, read2, trunk, expect_trunk, f"{outputFilePrefix}.__TRUNK__.2.fastq.gz", fill_length, compresslevel)

    trunk += expect_trunk

  logger.info("done")

def split_dynamic_single_end(logger, inputFiles, outputFilePrefix, min_file_size_gb, trunk_file_size_gb, fill_length, compresslevel=1, call_fastqsplitter=False): 
  trunk = 1
  for read1 in inputFiles:
    read1_gb = get_file_size(read1)
    file1 = f"{outputFilePrefix}.{str(trunk).zfill(fill_length)}.fastq.gz"

    if read1_gb < min_file_size_gb:
      logger.info(f"softlink {read1} => {file1}")
      do_symlink(read1, file1)
      trunk += 1
      continue

    expect_trunk = math.ceil(read1_gb / trunk_file_size_gb)

    if not call_fastqsplitter:
      split_by_trunk(logger, read1, False, outputFilePrefix, expect_trunk, trunk, fill_length, compresslevel)
    else:
      do_fastqsplitter(logger, read1, trunk, expect_trunk, f"{outputFilePrefix}.__TRUNK__.fastq.gz", fill_length, compresslevel)

    trunk += expect_trunk

  logger.info("done")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Scatter pairend big FASTQ files to small FASTQ files",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  #input can be multiple FASTQ files, such as F1_R1,F1_R2,F2_R1,F2_R2,F3_R1,F3_R2 in case some of them are more larger than others.
  #we don't know how many reads per sample
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', default="-", help="Output file prefix", required=NOT_DEBUG)
  parser.add_argument('--is_single_end', action='store', nargs='?', help="Is single end?")
  parser.add_argument('--min_file_size_gb', action='store', nargs='?', type=float, default=10, help="Min file size for split in GB")
  parser.add_argument('--trunk_file_size_gb', action='store', nargs='?', type=float, default=5, help="Expect file size for splitted file in GB")
  parser.add_argument('--fill_length', action='store', nargs='?', type=int, default=3, help="Trunk name length (fill with zero)")
  parser.add_argument('--compresslevel', action='store', nargs='?', type=int, default=1, help="Compress level, 1: fastest, 9: slowest")
  parser.add_argument('--call_fastqsplitter', action='store_true', help="Call fastqsplitter to speed up")
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-1_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-1_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-2_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-2_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-3_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-3_2.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-4_1.fq.gz,/scratch/cqs/breast_cancer_spore/clean/02-220/V350088726_L03_B5GHUMexlkRABOA-4_2.fq.gz"
    args.outputPrefix = "/scratch/cqs/shengq2/temp/batch02/bwa_00_splitFastq/result/P02-220"
    args.is_single_end = False
    args.min_file_size_gb = 10
    args.trunk_file_size_gb = 5
  
  logger = logging.getLogger('splitFastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  inputFiles = args.input.split(",")
  
  if args.is_single_end:
    split_dynamic_single_end(logger, inputFiles, args.outputPrefix, args.min_file_size_gb, args.trunk_file_size_gb, args.fill_length, args.compresslevel, args.call_fastqsplitter)
  else:
    split_dynamic_paired_end(logger, inputFiles, args.outputPrefix, args.min_file_size_gb, args.trunk_file_size_gb, args.fill_length, args.compresslevel, args.call_fastqsplitter)
  
if __name__ == "__main__":
    main()
