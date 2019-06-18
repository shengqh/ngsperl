#!/usr/bin/env python
# encoding: utf-8
"""
bravo_AF_filter.py
filter bravo_AF info.
"""
import argparse
import sys
import logging
import os
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

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="filter VCF",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NotDEBUG)
parser.add_argument('-t', '--threshold', action='store', default=0.001, type=int, nargs='?', help='AF value')
parser.add_argument('--debug', action='store_true', help="Output debug information", default=False)

args = parser.parse_args()

if DEBUG:
  args.input = "./IPF_batch2.smoove.square.anno.vcf.gz"
  args.output = "./IPF_batch2.smoove.square.anno.vcf.Filtered.vcf"


t_AF=float(args.threshold)


logger = initialize_logger(args.output + ".log", 'bravo_AF_filter', args.debug)
logger.info(str(args))
with open(args.output, "w") as fout:
    if args.input.endswith(".gz"):
      fin = gzip.open(args.input, 'rb')
    else:
      fin = open(args.input, "r")
    try:
      while True:
        line = fin.readline()
        if line.startswith("Chr"):
          fout.write(line)
          vcfheaders = line.rstrip().split("\t")
          info_index = vcfheaders.index("INFO")
          
          break
        else:
          fout.write(line)

      for line in fin:
        if "bravo_AF" in line:
            snv = line.rstrip().split('\t')
            info_parts = snv[info_index].split(";")
            bravo_af=info_parts[-1].split("=")[1]
            if float(bravo_af) <= t_AF:
                fout.write(line)
        else:
            fout.write(line)
          
      logger.info("Done")
    finally:
      fin.close()
