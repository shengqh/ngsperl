import argparse
import sys
import logging
import os
import os.path
import csv
import gzip
import re

def is_version_2():
  return sys.version_info[0] < 3

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

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Filter vcf samples.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output MAF file name", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_03_filterMAF/result/Ciombor_ExomeSeq.maf_filtered.vcf.discard"
  args.output = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_03_filterMAF/result/Ciombor_ExomeSeq.maf_filtered.vcf.discard.slim"

logger = initialize_logger(args.output + ".log", 'filterVcfSample', True)
logger.info(str(args))

outputTemp = args.output + ".tmp"
with open(outputTemp, "w") as fout:
  if args.input.endswith(".gz"):
    if is_version_2():
      fin = gzip.open(args.input, 'rb')
    else:
      fin = gzip.open(args.input, 'rt')
  else:
    fin = open(args.input, "r")
  try:
    sample_index = 9
    for line in fin:
      if line.find("#CHROM") != -1:
        continue    

      snv = line.rstrip().split('\t')
      fout.write("%s" % ("\t").join(snv[0:9]))
      
      for si in range(sample_index, len(snv)):
        sampleData = snv[si]
        if sampleData.startswith("./.:") or sampleData.startswith("0/0:") or sampleData.startswith("0|0"):
          continue
        fout.write("\t%s" % sampleData)
      fout.write("\n")
  finally:
    fin.close()
      
if os.path.isfile(args.output):
  os.remove(args.output)
os.rename(outputTemp, args.output)

logger.info("done.")