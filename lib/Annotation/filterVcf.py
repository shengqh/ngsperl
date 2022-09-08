import argparse
import sys
import logging
import os
import os.path
import csv
import gzip
import re
from Bio import SeqIO, bgzf

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

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Filter vcf by allele frequency. If SNV is marked as 0/1 and more than assigned percentage (default 90%) samples with minor allele frequency less than assigned MAF (default 0.3), the SNV will be discarded. This criteria applies to both high quality samples (with min_depth, min_genotype_quality) first, then applies to all samples.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-f', '--frequency', action='store', default=0.3, type=float, nargs='?', help='Min minor allele frequency for valid genotype 1')
parser.add_argument('-p', '--percentage', action='store', default=0.9, type=float, nargs='?', help='Max sample percentage allowed with MAF less than frequency defined')
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output MAF file name", required=NotDEBUG)
parser.add_argument('--min_inbreeding_coeff', action='store', default=-0.2, type=float, nargs='?', help='Min inbreeding coefficient')
parser.add_argument('--min_depth', action='store', default=10, type=int, nargs='?', help='Min depth in at least one sample')
parser.add_argument('--min_genotype_quality', action='store', default=20, type=int, nargs='?', help='Min quality in at least one sample')
parser.add_argument('--SNP_only', action='store_true', help='Keep SNP only, remove indel SNV')
parser.add_argument('--debug', action='store_true', help="Output debug information", default=False)

args = parser.parse_args()
if DEBUG:
  args.input = "/gpfs23/scratch/cqs/fanr1/202111_ICARE_BEST/DNA/20211207_DNAseq_Analysis/bwa_g4_refine_gatk4_SNV_02_vqsr/result/human_exomeseq.indels.snp.recal.pass.norm.nospan.vcf.gz"
  args.output = "/scratch/cqs/shengq2/temp/human_exomeseq.maf_filtered.vcf.gz"


percentage=float(args.percentage)
frequency=float(args.frequency)

logger = initialize_logger(args.output + ".log", 'filterVcf', args.debug)
logger.info(str(args))

basename = os.path.splitext(args.output)[0]

if args.output.endswith(".gz"):
  outputTemp = basename + ".tmp.gz"
  fout = bgzf.BgzfWriter(outputTemp, "wb")
  fdiscard = bgzf.BgzfWriter(basename + ".discard.gz", "wb")
else:
  outputTemp = basename + ".tmp"
  fout = open(outputTemp, "wt")
  fdiscard = open(basename + ".discard", "wt")

if args.input.endswith(".gz"):
  if is_version_2():
    fin = gzip.open(args.input, 'rb')
  else:
    fin = gzip.open(args.input, 'rt')
else:
  fin = open(args.input, "r")

with fout:
  with fdiscard:
    with fin:
      while True:
        line = fin.readline()
        if line.find("#CHROM") != -1:
          fout.write(line)
          fdiscard.write(line)
          vcfheaders = line.rstrip().split("\t")
          format_index = vcfheaders.index("FORMAT")
          info_index = vcfheaders.index("INFO")
          sample_index = format_index + 1
          AD_index = -1
          GQ_index = -1
          DP_index = -1
          break
        else:
          fout.write(line)
      
      totalsnv = 0
      FailInbreedingCoeff = 0
      FailQuality = 0
      FailPercentage = 0
      indel = 0
      lastChrom = ""

      for line in fin:
        snv = line.rstrip().split('\t')
        
        totalsnv = totalsnv + 1
        if totalsnv % 10000 == 0:
          logger.info("Chrom=%s,Processed=%d,FailInbreedingCoeff=%d,FailQuality=%d,FailPercentage=%d,Passed=%d" % (snv[0], totalsnv, FailInbreedingCoeff, FailQuality, FailPercentage, totalsnv - FailInbreedingCoeff - FailQuality - FailPercentage))
        
        if args.SNP_only and (len(snv[3]) > 1 or len(snv[4]) > 1):
          fdiscard.write("Indel=%s" % (line))
          indel += 1
          continue

        gt1count = 0
        gt1lessCount = 0
        highQuality_gt1count = 0
        highQuality_gt1lessCount = 0
        if AD_index == -1:
          format_parts = snv[format_index].split(":")
          AD_index = format_parts.index("AD")
          GQ_index = format_parts.index("GQ")
          DP_index = format_parts.index("DP")

        info = snv[info_index]
        m = re.search('InbreedingCoeff=(.+?);', info)
        if m:
          ic = float(m.group(1))
          if ic < args.min_inbreeding_coeff:
            fdiscard.write("FailInbreedingCoeff=%f : %s" % (ic, line))
            FailInbreedingCoeff = FailInbreedingCoeff + 1
            continue
        
        highQuality = False
        highQualityGt1 = list()
        for si in range(sample_index, len(vcfheaders)):
          sampleData = snv[si]
          if sampleData == './.' or sampleData.startswith("0/0:") or sampleData.startswith("0|0"):
            continue
          
          is_current_highQuality = False
          parts = sampleData.split(":")
          if len(parts) <= GQ_index or len(parts) <= DP_index:
            logger.error('failed sample data %s in line %s' % (sampleData, line))
          if parts[GQ_index] != '.' and parts[DP_index] != '.':
            gq = int(parts[GQ_index])
            dp = int(parts[DP_index])
            if gq >= args.min_genotype_quality and dp >= args.min_depth:
              is_current_highQuality = True
              highQuality = True

          is_gt1 = sampleData.startswith("0/")
          if is_gt1:
            split_gt1 = "/"
          else:
            is_gt1 = sampleData.startswith("0|")
            if is_gt1:
              split_gt1 = "|"

          if is_gt1:
            gt1count += 1
            if is_current_highQuality:
              highQuality_gt1count += 1
            
            gts = parts[0].split(split_gt1)
            gt = int(gts[1])
            
            adList = parts[AD_index].split(",")
            ad0 = float(adList[0])
            ad1 = float(adList[gt])
            
            if ad0 + ad1 == 0:
              gt1lessCount = gt1lessCount + 1
              if is_current_highQuality:
                highQuality_gt1lessCount += 1
            elif(ad1  / (ad0 + ad1) < frequency):
              gt1lessCount = gt1lessCount + 1
              if is_current_highQuality:
                highQuality_gt1lessCount += 1
        
        if not highQuality:
          FailQuality = FailQuality + 1
          fdiscard.write("lowQuality\t%s\n" %(line.rstrip()))
          continue 

        if gt1count == 0:
          logger.debug("No GT1:" + line)
          fout.write(line)
          continue

        #check high quality GT1 
        if highQuality_gt1count > 0:          
          gt1lessPercentage = highQuality_gt1lessCount * 1.0 / highQuality_gt1count
          if gt1lessPercentage >= percentage:
            FailPercentage = FailPercentage + 1
            fdiscard.write("failedCount/gt1count of highQuality GT1=%d/%d\t%s\n" %(highQuality_gt1lessCount, highQuality_gt1count, line.rstrip()))
            continue

        #check all GT1           
        gt1lessPercentage = gt1lessCount * 1.0 / gt1count
        if gt1lessPercentage < percentage:
          fout.write(line)
        else:
          FailPercentage = FailPercentage + 1
          fdiscard.write("failedCount/gt1count of all GT1=%d/%d\t%s\n" %(gt1lessCount, gt1count, line.rstrip()))
      
      logger.info("Total\t%d" % (totalsnv))
      logger.info("Indel\t%d" % (indel))
      logger.info("FailInbreedingCoeff\t%d" % (FailInbreedingCoeff))
      logger.info("FailQuality\t%d" % (FailQuality))
      logger.info("FailPercentage\t%d" % (FailPercentage))
      logger.info("Done")

      with open(args.output + ".summary.txt", "wt") as fs:
        fs.write("Category\tCount\n")
        fs.write("Total\t%d\n" % (totalsnv))
        fs.write("Indel\t%d\n" % (indel))
        fs.write("FailInbreedingCoeff\t%d\n" % (FailInbreedingCoeff))
        fs.write("FailQuality\t%d\n" % (FailQuality))
        fs.write("FailPercentage\t%d\n" % (FailPercentage))
        fs.write("Pass\t%d\n" % (totalsnv - indel - FailInbreedingCoeff - FailQuality - FailPercentage))
      
if os.path.isfile(args.output):
  os.remove(args.output)
os.rename(outputTemp, args.output)

if args.output.endswith(".gz"):
  cmd = "tabix %s" % (args.output)
  logger.info(cmd)
  os.system(cmd)

logger.info("done.")
