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

def isSexChromosome(chrom):
  return chrom == "X" or chrom == "Y" or chrom == "chrX" or chrom == "chrY"

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

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Calculate IBS",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-f', '--family', action='store', nargs='?', help="Input family info file (two columns, first is sample, second is family, without header)", required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output IBS score name", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_08_IBS/result/test.vcf"
  args.output = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_08_IBS/result/test.py.score.txt"
  args.family = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_08_IBS/result/Ciombor_ExomeSeq__fileList1.list"
  # args.input = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_03_filterMAF/result/Ciombor_ExomeSeq.maf_filtered.vcf"
  # args.output = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_08_IBS/result/Ciombor_ExomeSeq.py.score.txt"
  # args.family = "/scratch/cqs/ramirema/20190610_Ciombior_ExomeSeq/results/bwa_refine_gatk4_SNV_08_IBS/result/Ciombor_ExomeSeq__fileList1.list"

logger = initialize_logger(args.output + ".log", 'IBS', True)
logger.info(str(args))

gtscore = {
  "00": {
    "00":-1,
    "01":1,
    "10":1,
    "11":0,
  },
  "01": {
    "00":1,
    "01":2,
    "10":2,
    "11":1,
  },
  "10": {
    "00":1,
    "01":2,
    "10":2,
    "11":1,
  },
  "11": {
    "00":0,
    "01":1,
    "10":1,
    "11":2,
  },
}

#print(gtscore)
vcfheaders = []
sample_index = 0
scoreMap = {}

if args.input.endswith(".gz"):
  if is_version_2():
    fin = gzip.open(args.input, 'rb')
  else:
    fin = gzip.open(args.input, 'rt')
else:
  fin = open(args.input, "r")
try:
  sample_index = 9
  while(True):
    line = fin.readline()
    if line.find("#CHROM") != -1:
      vcfheaders = line.rstrip().split("\t")
      h = len(vcfheaders)
      format_index = vcfheaders.index("FORMAT")
      info_index = vcfheaders.index("INFO")
      sample_index = format_index + 1
      sample_names = vcfheaders[sample_index:]
      GQ_index = -1
      DP_index = -1
      MAX_index = -1
      all_range = range(sample_index, len(vcfheaders))
      s1_range = range(sample_index, len(vcfheaders)-1)
      for s1 in s1_range:
        s1map = {}
        scoreMap[s1] = s1map
        for s2 in range(s1+1, h):
          s1map[s2] = [0,0,0]
      break    
  
  snvcount = 0
  for line in fin:
    snvcount += 1
    if snvcount % 1000 == 0:
      logger.info("Processed %d ..." % snvcount)

    #if snvcount % 500 == 0:
      #print(scoreMap)
    #  break

    snv = line.rstrip().split('\t')
    chrom = snv[0]
    if isSexChromosome(chrom):
      continue
    
    if GQ_index == -1:
      format_parts = snv[format_index].split(":")
      GQ_index = format_parts.index("GQ")
      DP_index = format_parts.index("DP")
      MAX_index = max(GQ_index, DP_index) + 1

    validIndecies = []
    for idx in all_range:
      sdata = snv[idx]
      if genotypeValid(sdata, MAX_index, GQ_index, DP_index, 21, 6):
        validIndecies.append(idx)
        snv[idx] = sdata[0] + sdata[2]

    #for idx in validIndecies:
    #  print("%d : %s" % (idx, snv[idx]))

    for v1 in range(0, len(validIndecies)-1):
      s1 = validIndecies[v1]
      s1gt = snv[s1]
      #if s1gt not in gtscore:
      #  raise Exception("Genotype %s not defined in gtscore" % s1gt)

      s1map = gtscore[s1gt]
      for v2 in range(v1+1, len(validIndecies)):
        s2 = validIndecies[v2]
        s2gt = snv[s2]
        # if s2gt not in s1map:
        #   raise Exception("Genotype %s not defined in map %s" % (s2gt, s1gt))

        score = s1map[s2gt]
        if score == -1:
          continue

        mdata = scoreMap[s1][s2]
        mdata[score] += 1
finally:
  fin.close()

basename = os.path.splitext(args.output)[0]

with open(args.output, "w") as fout:
  fout.write(",%s\n" % ",".join(vcfheaders[sample_index:]))
  for fi in all_range:
    fout.write(vcfheaders[fi])
    for fj in all_range:
      if (fj <= fi):
        fout.write(",NA")
      else:
        mdata = scoreMap[fi][fj]
        fout.write(",%d:%d:%d" % (mdata[0], mdata[1], mdata[2] ) )
    fout.write("\n")
      
mean_file = basename+ ".mean.csv"
with open(mean_file, "w") as fout:
  fout.write(",%s\n" % ",".join(vcfheaders[sample_index:]))
  for fi in all_range:
    fout.write(vcfheaders[fi])
    for fj in all_range:
      if (fj <= fi):
        fout.write(",NA")
      else:
        mdata = scoreMap[fi][fj]
        total = mdata[0] + mdata[1] + mdata[2]
        if total == 0:
          fout.write(",NA")
        else:
          meanvalue = (mdata[1] + mdata[2] * 2) / (mdata[0] + mdata[1] + mdata[2])
          fout.write(",%.3lf" % meanvalue )
    fout.write("\n")

  realpath = os.path.dirname(os.path.realpath(__file__))
  rPath = realpath + "/IBS.R"

  cmd = "Rscript %s -i %s -f %s -o %s" % (rPath, mean_file, args.family, basename)
  logger.info(cmd)
  os.system(cmd)

logger.info("done.")