#!/usr/bin/env python3

import argparse
import sys
import logging
import os
import csv
import gzip
  
DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Convert TOPMed database to annovar format.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input TOPMed file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output Output file", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/references/TOPMed/bravo-dbsnp-all.vcf.gz"
  args.output = "/scratch/cqs/references/annovar/humandb/hg38_topmed05.txt"

logger = logging.getLogger('topmed2annovar')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

if args.input.endswith(".gz"):
  fin = gzip.open(args.input, 'rt')
else:
  fin = open(args.input, "rt")
  
with fin:
  tmpFile = args.output + ".tmp"
  with open(tmpFile, "w") as fout:
    fout.write("#Chr\tStart\tEnd\tRef\tAlt\tTOPMed\n")
    for line in fin:
      if line.startswith('#'):
        continue
      parts = line.rstrip().split("\t")
      chr = parts[0]
      if chr.startswith("chr"):
        chr = chr[3:]
      start = parts[1]
      end = parts[1]
      ref = parts[3]
      alt = parts[4]
      info = parts[7]
      infoparts = info.split(';')
      af = ""
      for ip in infoparts:
        ipparts = ip.split("=")
        if ipparts[0] == "AF":
          af = ipparts[1]
          break
      if af == "":
        continue
      fout.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chr, start, end, ref, alt, af))
    
  if os.path.isfile(args.output):
    os.remove(args.output)
  os.rename(tmpFile, args.output)

  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/annovar_idx.pl"
  cmd = "perl " + rPath + " " + args.output + " 1000 > " + args.output + ".idx"
  logger.info(cmd)
  os.system(cmd)
          
  logger.info("done.")

