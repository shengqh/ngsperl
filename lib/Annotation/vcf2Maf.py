import argparse
import sys
import logging
import os
import csv
from IPython.utils import ulinecache

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Post vcf to maf processor.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input MAF file', required=NotDEBUG)
parser.add_argument('-v', '--vcf', action='store', nargs='?', help='Input VCF file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output MAF file name", required=NotDEBUG)

args = parser.parse_args()
if DEBUG:
  args.input = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_hc_gvcf_hardfilter_vep/result/linton_exomeseq_2118.pass.vcf.maf.tmp.tmp"
  args.vcf = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_hc_gvcf_hardfilter/result/linton_exomeseq_2118.pass.vcf.tmp"
  args.out = "T:/Shared/Labs/Linton Lab/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_hc_gvcf_hardfilter_vep/result/linton_exomeseq_2118.pass.vcf.tmp.maf.txt"

logger = logging.getLogger('vcf2maf')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

with open(args.vcf, "r") as fin:
  while True:
    line = fin.readline()
    if "#CHROM" in line:
      vcfheaders = line.rstrip().split("\t")
      sample_index = vcfheaders.index("FORMAT") + 1
      break
  vcfdata = [row.split('\t') for row in fin]

with open(args.input, "r") as fin:
  version = fin.readline()
  mafheaders = fin.readline().split("\t")
  Tumor_Sample_Barcode_index = mafheaders.index("Tumor_Sample_Barcode")
  t_depth_index = mafheaders.index("t_depth")
  t_ref_count_index = mafheaders.index("t_ref_count")
  t_alt_count_index = mafheaders.index("t_alt_count")
  mafdata = [row.split('\t') for row in fin]

if len(vcfdata) != len(mafdata):
  print("ERROR, there are %d data in vcf while %d data in maf\n" % (len(vcfdata), len(mafdata)))
else:
  with open(args.out, "w") as fout:
    fout.write(version)
    fout.write('\t'.join(mafheaders))
    for idx, mafitem in enumerate(mafdata):
      vcfitem = vcfdata[idx]
      for sample_idx in range(sample_index, len(vcfheaders)):
        sample_name = vcfheaders[sample_idx]
        sample_data = vcfitem[sample_idx]
        if "./.:" in sample_data:
          continue
        if "0/0:" in sample_data:
          continue
        parts = sample_data.split(":")
        alleles = parts[1].split(",")
        mafitem[Tumor_Sample_Barcode_index] = sample_name
        mafitem[t_depth_index] = parts[2]
        mafitem[t_ref_count_index] = alleles[0]
        mafitem[t_alt_count_index] = alleles[1]
        fout.write("\t".join(mafitem))
