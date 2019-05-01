#!/usr/bin/env python
# encoding: utf-8
"""
SV_filter.py
filter SVs in vcf format with read depth in family-wise.
"""
import argparse
import sys
import logging
import os
import csv
import gzip
import re
from peds import open_ped, Family, Person

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

parser = argparse.ArgumentParser(description="filter SV that produced by smoove",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output file name", required=NotDEBUG)
parser.add_argument('-p', '--pedigree', action='store', nargs='?', help="Input pedigree file")
parser.add_argument('--max_DHFFC_DEL', action='store', default=0.7, type=int, nargs='?', help='Max DHFFC value for deletion in at least one sample')
parser.add_argument('--min_DHFFC_DUP', action='store', default=1.25, type=int, nargs='?', help='Min DHFFC value for duplication in at least one sample')
parser.add_argument('--debug', action='store_true', help="Output debug information", default=False)

args = parser.parse_args()

if DEBUG:
  args.input = "./IPF_batch2.smoove.square.anno.vcf.gz"
  args.output = "./IPF_batch2.smoove.square.anno.vcf.Filtered.vcf"


t_del=float(args.max_DHFFC_DEL)
t_dup=float(args.min_DHFFC_DUP)


logger = initialize_logger(args.output + ".log", 'SV_filter', args.debug)
logger.info(str(args))


ped={}
family={}

if args.pedigree:
    logger.info("You provied the Family file.")
    with open(args.pedigree) as handle:
        sep = "\t"
        for line in handle:
            # ignore header and comment lines
            if line.startswith('family_id{}'.format(sep)) or line.startswith('#'):
                continue
            record=line.split(sep, 2)
            fam_id=record[0]
            person_id=record[1]
            ped[person_id]=fam_id
            if fam_id not in family:
              family[fam_id] = 1
            else:
              family[fam_id] += 1

with open(args.output, "w") as fout:
    if args.input.endswith(".gz"):
      fin = gzip.open(args.input, 'rb')
    else:
      fin = open(args.input, "r")
    try:
      while True:
        line = fin.readline()
        if "#CHROM" in line:
          fout.write(line.rstrip())
          vcfheaders = line.rstrip().split("\t")
          filter_index = vcfheaders.index("FILTER")
          format_index = vcfheaders.index("FORMAT")
          info_index = vcfheaders.index("INFO")
          sample_index = format_index + 1
          
          break
        else:
          fout.write(line)
      fout.write("\t"+"SummuryByFamily"+"\t"+"NumAnySample"+"\t"+"NumAnyFamily"+"\t"+"NumAllFamily"+"\n")
      totalsnv = 0
      for line in fin:
        score={}
        snv = line.rstrip().split('\t')
        totalsnv = totalsnv + 1
        highQuality = False
        info_parts = snv[info_index].split(";")
        SV_type=info_parts[0].split("=")[1]
        if info_parts[-1].startswith("smoove_gene"):
            SV_MSHQ=info_parts[-2].split("=")[1]
        else:
            SV_MSHQ=info_parts[-1].split("=")[1]
#        print SV_MSHQ
        if SV_MSHQ <4:
          continue
        if SV_type != "INV" and SV_type != "BND":
          if SV_type == "DEL":
            format_parts = snv[format_index].split(":")
            GT_index = format_parts.index("GT")
            DHFFC_index = format_parts.index("DHFFC")
            #SHQ_index = format_parts.index("SHQ")
            nIndex = format_index
            for si in range(sample_index, len(vcfheaders)):
              nIndex += 1
              sampleData = snv[si]
              if sampleData.startswith("0/0:") or sampleData.startswith("0|0"):
                continue
              parts = sampleData.split(":")
              if parts[GT_index].startswith("0/") or parts[GT_index].startswith("0|") :
                shq = int(parts[-1])
                dhffc = float(parts[DHFFC_index])
                if dhffc <= t_del and shq == 4:
                    highQuality = True
                    sampleName=vcfheaders[nIndex]
                    print  sampleName, ped[sampleName]
                    if ped[sampleName] not in score:
                        score[ped[sampleName]]=1
                    else:
                        score[ped[sampleName]] += 1

          elif SV_type == "DUP":
            format_parts = snv[format_index].split(":")
            GT_index = format_parts.index("GT")
            DHFFC_index = format_parts.index("DHFFC")
            #SHQ_index = format_parts.index("SHQ")
            nIndex=format_index
            for si in range(sample_index, len(vcfheaders)):
              nIndex += 1
              sampleData = snv[si]
              if sampleData.startswith("0/0:") or sampleData.startswith("0|0"):
                continue
              parts = sampleData.split(":")
              if parts[GT_index].startswith("0/") or parts[GT_index].startswith("0|") :
                shq = int(parts[-1])
                dhffc = float(parts[DHFFC_index])
                if dhffc >= t_dup and shq == 4:
                    highQuality = True
                    sampleName=vcfheaders[nIndex]
                    print sampleName
                    if ped[sampleName] not in score:
                        score[ped[sampleName]]=1
                    else:
                        score[ped[sampleName]] += 1

          if highQuality == True:
            snv[filter_index]="PASS"
            frontline='\t'.join(map(str, snv))
            bline=[]
            total=0
            count=0
            for (key,val) in score.iteritems():
                s_pec=val*1.0/family[key]
                count += val
                if val == family[key]:
                    total += 1
                bline.append('%s=%.2f' % (key,s_pec))
            backline=':'.join(bline)
            fout.write(frontline+"\t"+backline+"\t"+str(count)+"\t"+str(len(score))+"\t"+str(total)+"\n")
        #  else:
        #    fout.write(line)
        #else:
        #  fout.write(line)
          
      logger.info("Done")
    finally:
      fin.close()