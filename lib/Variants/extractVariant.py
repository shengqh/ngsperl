import argparse
import sys
import logging
import os
import os.path
import csv
import gzip
import re

def initialize_logger(logfile, logname, isDebug):
  logger = logging.getLogger(logname)
  loglevel = logging.DEBUG if isDebug else logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

def extract_variant(logger, fout, vcf_name, vcf_file, snp_chrom, snp_position):
  logger.info(f"Processing {vcf_name} ...")
  if not os.path.exists(vcf_file):
    return
  
  if vcf_file.endswith(".gz"):
    fin = gzip.open(vcf_file, 'rt')
  else:
    fin = open(vcf_file, "r")

  with fin:
    sample_index = 9
    while(True):
      line = fin.readline()
      if line.startswith("#CHROM"):
        break
    
    snvcount = 0
    for line in fin:
      snvcount += 1
      if snvcount % 1000000 == 0:
        logger.info("Processed %d ..." % snvcount)

      snv = line.rstrip().split('\t', 2)
      chrom = snv[0]
      if chrom != snp_chrom:
        continue
      
      position = snv[1]
      if position != snp_position:
        continue
      
      fout.write(f"{vcf_name}\t{line}")
      break

def extract_variant_list(logger, listFile, outputFile, snp_chrom, snp_position):
  with open(listFile, "rt") as fin, open(outputFile, "wt") as fout:
    fout.write("Sample\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGenotype\n")
    for line in fin:
      parts = line.rstrip().split('\t')
      vcf_name = parts[1]
      vcf_file = parts[0]
      extract_variant(logger, fout, vcf_name, vcf_file, snp_chrom, snp_position)

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Calculate Chromosome Count",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF list file', required=NotDEBUG)
parser.add_argument('--snp_chrom', action='store', help="Input snp chrom", required=NotDEBUG)
parser.add_argument('--snp_position', action='store', help="Input snp position", required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NotDEBUG)

args = parser.parse_args()

logger = initialize_logger(args.output + ".log", 'extract_snp', True)
logger.info(str(args))

extract_variant_list(logger, args.input, args.output, args.snp_chrom, args.snp_position)

logger.info("done.")