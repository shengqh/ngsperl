import argparse
import sys
import logging
import os
import gzip
import math
from _ctypes import ArgumentError

def should_remove(sequence, removeSequences):
  for seq in removeSequences:
    if seq in sequence:
      return(True)
  return(False)

def remove(logger, inputFiles, is_single_end, outputFilePrefix, removeSequences): 
  if is_single_end:
    read1files = inputFiles
    read1fout = gzip.open(outputFilePrefix + ".fastq.gz", "wt")
  else:
    read1files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 0]
    read2files = [inputFiles[idx] for idx in range(0, len(inputFiles)) if idx % 2 == 1]
    read1fout = gzip.open(outputFilePrefix + ".1.fastq.gz", "wt")
    read2fout = gzip.open(outputFilePrefix + ".2.fastq.gz", "wt")

  for idx in range(0, len(read1files)):
    logger.info("reading " + read1files[idx] + " ...")
    read1fin = gzip.open(read1files[idx], "rt")
    if not is_single_end:
      logger.info("reading " + read2files[idx] + " ...")
      read2fin = gzip.open(read2files[idx], "rt")

    read_count = 0
    removed_read1_count = 0
    removed_read2_count = 0
    while(True):
      r1line1 = read1fin.readline()
      if not r1line1:
        break

      if r1line1[0] != '@':
        raise Exception("Wrong id: %s" % r1line1)
      
      r1line2 = read1fin.readline()
      r1line3 = read1fin.readline()
      r1line4 = read1fin.readline()

      removed = should_remove(r1line2, removeSequences)
      if removed:
        removed_read1_count += 1

      if not is_single_end:
        r2line1 = read2fin.readline()
        if not r2line1:
          break

        if r2line1[0] != '@':
          raise Exception("Wrong id: %s" % r2line1)
        
        r2line2 = read2fin.readline()
        r2line3 = read2fin.readline()
        r2line4 = read2fin.readline()

        if not removed:
          removed = should_remove(r2line2, removeSequences)
          if removed:
            removed_read2_count += 1


      read_count += 1
      if read_count % 100000 == 0:
        total_removed = removed_read1_count + removed_read2_count
        logger.info(f"processed {read_count} reads, removed {total_removed} ...")
    
      if removed:
        continue
    
      read1fout.write(r1line1)
      read1fout.write(r1line2)
      read1fout.write(r1line3)
      read1fout.write(r1line4)

      if not is_single_end:
        read2fout.write(r2line1)
        read2fout.write(r2line2)
        read2fout.write(r2line3)
        read2fout.write(r2line4)

    read1fin.close()
    if not is_single_end:
      read2fin.close()

  read1fout.close()
  if not is_single_end:
    read2fout.close()

  with open(outputFilePrefix + ".info", "wt") as fout:
    fout.write("read\tcount\n")
    fout.write(f"total\t{read_count}\n")
    fout.write(f"remove_read1\t{removed_read1_count}\n")
    if not is_single_end:
      fout.write(f"remove_read2\t{removed_read2_count}\n")

  logger.info("done")

def main():
  DEBUG = False
  NOT_DEBUG = not DEBUG
  
  parser = argparse.ArgumentParser(description="Remove reads with specific sequences",
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  
  parser.add_argument('-i', '--input', action='store', nargs='?', help='Input Fastq files (first,second for pairend data)', required=NOT_DEBUG)
  parser.add_argument('-o', '--outputPrefix', action='store', nargs='?', help="Output file prefix", required=NOT_DEBUG)
  parser.add_argument('-s', '--sequence', action='store', nargs='?', help="Input specific sequences (split by ',')", required=NOT_DEBUG)
  parser.add_argument('--is_single_end', action='store', nargs='?', help="Is single end?")
  
  args = parser.parse_args()
  
  if DEBUG:
    args.input = "/scratch/jbrown_lab/shengq2/projects/20201210_rnaseq_PRJNA627465_mm10_LMNA/sra2fastq/result/Untreated_Parent167_3_1.fastq.gz,/scratch/jbrown_lab/shengq2/projects/20201210_rnaseq_PRJNA627465_mm10_LMNA/sra2fastq/result/Untreated_Parent167_3_2.fastq.gz"
    args.is_single_end = False
    args.outputPrefix = "/scratch/jbrown_lab/shengq2/projects/20201210_rnaseq_PRJNA627465_mm10_LMNA/remove_PCR_primer/Untreated_Parent167_3_1"
    args.sequence = "GTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCC"
  
  logger = logging.getLogger('removeSequence')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  
  inputFiles = args.input.split(",")
  removeSequences = args.sequence.split(",")
  
  remove(logger, inputFiles, args.is_single_end, args.outputPrefix, removeSequences)
  
if __name__ == "__main__":
    main()
