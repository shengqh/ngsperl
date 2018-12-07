import argparse
import sys
import logging
import os
import csv

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Filter vcf by allele frequency. If SNV is marked as 0/1 and more than assigned percentage (default 90%) samples with minor allele frequency less than assigned MAF (default 0.3), the SNV will be discarded.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-p', '--percentage', action='store', default=0.9, nargs='?', help='Max sample percentage allowed')
parser.add_argument('-f', '--frequency', action='store', default=0.3, nargs='?', help='Max minor allele frequency')
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output MAF file name", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_hc_gvcf_hardfilter/result/linton_exomeseq_2118.pass.vcf"
  args.output = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_hc_gvcf_hardfilter_vep/result/linton_exomeseq_2118.pass.filtered.vcf"

logger = logging.getLogger('filterVcf')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')


with open(args.output, "w") as fout:
  with open(args.output + ".discard", "w") as fdiscard:
    with open(args.input, "r") as fin:
      while True:
        line = fin.readline()
        if "#CHROM" in line:
          fout.write(line)
          vcfheaders = line.rstrip().split("\t")
          format_index = vcfheaders.index("FORMAT")
          sample_index = format_index + 1
          AD_index = -1
          break
        else:
          fout.write(line)
      
      totalsnv = 0
      svnzero = 0
      for line in fin:
        snv = line.split('\t')
        
        totalsnv = totalsnv + 1
        
        gt1count = 0
        gt1lessCount = 0
        if AD_index == -1:
          AD_index = snv[format_index].split(":").index("AD")
        
        haszero = False
        for si in range(sample_index, len(vcfheaders)):
          sampleData = snv[si]
          if "0/0:" in sampleData:
            continue
          
          if "0/" in sampleData:
            gt1count = gt1count + 1
            
            parts = sampleData.split(":")
            gts = parts[0].split("/")
            gt = int(gts[1])
            
            adList = parts[AD_index].split(",")
            ad0 = float(adList[0])
            ad1 = float(adList[gt])
            
            if ad0 + ad1 == 0:
              haszero = True
              gt1lessCount = gt1lessCount + 1
            elif(ad1  / (ad0 + ad1) < args.frequency):
              gt1lessCount = gt1lessCount + 1
        
        if haszero:
          print("haszero: gt1count %d ~ failedCount %d : %s" %(gt1count, gt1lessCount, line.rstrip()))
          svnzero = svnzero + 1
          
        if gt1count == 0:
          fout.write(line)
          continue
          
        gt1lessPercentage = gt1lessCount * 1.0 / gt1count
        if gt1lessPercentage < args.percentage:
          fout.write(line)
        else:
          fdiscard.write("gt1count %d ~ failedCount %d\t%s\n" %(gt1count, gt1lessCount, line.rstrip()))
      print("haszero: %d out of %d" % (svnzero, totalsnv))
