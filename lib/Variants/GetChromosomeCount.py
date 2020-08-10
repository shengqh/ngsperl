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

def genotypeValid(sdata, MAX_index, GQ_index, DP_index, min_gq, min_dp):
  if sdata.startswith("."):
    return(False)

  parts = sdata.split(':', MAX_index)
  gq = int(parts[GQ_index])
  if (gq < min_gq):
    return(False)

  dp = int(parts[DP_index])
  if (dp < min_dp):
    return(False)

  return(True)

def getChromsomeCount(logger, vcfFile, outputFile=None):
  result = {}
  if vcfFile.endswith(".gz"):
    if is_version_2():
      fin = gzip.open(vcfFile, 'rb')
    else:
      fin = gzip.open(vcfFile, 'rt')
  else:
    fin = open(vcfFile, "r")
  try:
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
      if chrom in result:
        result[chrom] += 1
      else:
        result[chrom] = 1        
  finally:
    fin.close()

  if outputFile != None:
    with open(outputFile, "wt") as fout:
      fout.write("Chromosome\t%s\n" % name)
      for chrom in sorted(result.keys()):
        fout.write("%s\t%d\n" % (chrom, result[chrom]))

  return(result)


def getChromsomeCountList(logger, listFile, outputFile=None):
  result = []
  with open(args.input, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      name = parts[1]
      filepath = parts[0]
      curres = getChromsomeCount(logger, filepath)
      result.append(curres)
  
  if outputFile != None:
    with open(outputFile, "wt") as fout:
      fout.write("Chromosome\t%s\n" % name)
      for chrom in sorted(result.keys()):
        fout.write("%s\t%d\n" % (chrom, result[chrom]))

      files.append([parts[0], parts[1]])

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Calculate Chromosome Count",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('--is_list', action='store_true', help="Input file is list of vcf file", required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output IBS score name", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  #args.input = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_01_hc_gvcf/result/Tumor_99-24902.g.vcf"
  #args.output = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_01_hc_gvcf/result/Tumor_99-24902.g.vcf.chromosome"
  args.input = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_02_vqsr_4_gather/result/Ciombor_ExomeSeq.indels.snp.recal.pass.norm.nospan.vcf.gz"
  args.name = "test"
  args.output = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_02_vqsr_4_gather/result/Ciombor_ExomeSeq.indels.snp.recal.pass.norm.nospan.vcf.gz.chromosome"

logger = initialize_logger(args.output + ".log", 'IBS', True)
logger.info(str(args))

if args.is_list:
  getChromsomeCountList(logger, args.input, args.output)
else:
  getChromsomeCount(logger, args.input, args.output)

logger.info("done.")