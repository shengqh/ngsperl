import argparse
import sys
import logging
import os
import csv
import gzip
from subprocess import call

#https://software.broadinstitute.org/gatk/documentation/article.php?id=6926

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Add END description in vcf file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NotDEBUG, help='Input VCF file')
parser.add_argument('-o', '--output', action='store', nargs='?', required=NotDEBUG, help='Output VCF file')

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/ramirema/other_projects/20201215_Kristy_WES/bwa_refine_muTect2indel/result/H11_21642N_P10_1406T.somatic.indel.pass.vcf.gz"
  args.output = "/scratch/cqs/shengq2/temp/temp.vcf.gz"

print(args)

logger = logging.getLogger('addEND')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

logger.info("add END to vcf ...")
if args.output.endswith(".gz"):
  vcf_file = args.output[:-3]
else:
  vcf_file = args.output

with open(vcf_file, "wb") as fout:
  if args.input.endswith(".gz"):
    fin = gzip.open(args.input, 'rb')
    key = b"##INFO="
  else:
    fin = open(args.input, "r")
    key = "##INFO="
  try:
    bFirst = True
    for line in fin:
      if bFirst and line.startswith(key):
        fout.write(b'##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n')
        bFirst = False
      fout.write(line)
  finally:
    fin.close()

if args.output.endswith(".gz"):
  os.system("bgzip -f " + vcf_file)
  os.system("tabix -p vcf " + args.output)
    
logger.info("done.")
