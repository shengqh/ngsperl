import argparse
import sys
import logging
import os
import subprocess
import re
from Bio import SeqIO
from shutil import copyfile

DEBUG=True
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Merge impute2 result files",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', help='Input impute2 list file', required=NotDEBUG)
parser.add_argument('-o', '--output', action='store', nargs='?', help="Output plink file", required=NotDEBUG)

args = parser.parse_args()

if DEBUG:
  #rootFolder = "H:/shengquanhu/projects/macrae_linton/20190411_linton_megachip_2118_human/"
  rootFolder = "/scratch/cqs/shengq2/macrae_linton/20190411_linton_megachip_2118_human/"
  args.input = rootFolder + "impute2_merge/result/linton_exomeseq_2118__fileList1.list"
  args.sample = rootFolder + "plinkqc/result/linton_clean_gen.sample"
  args.output = rootFolder + "impute2_merge/result/linton_exomeseq_2118"

logger = logging.getLogger('mergeImpute2')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

class SNP:
  def __init__(self, chrom, name, position, firstAllele, secondAllele, genoTypes):
    self.chrom = chrom
    self.name = name
    self.position = position
    self.firstAllele = firstAllele
    self.secondAllele = secondAllele
    self.genoTypes = genoTypes

def readImpute2File(inputFile, chrom):
  result = []
  with open(inputFile, "r") as fin:
    for line in fin:
      parts = line.rstrip().split(' ')
      name = parts[1]
      position = int(parts[2])
      firstAllele = parts[3]
      secondAllele = parts[4]
      genoTypes = parts[5:]
      result.append(SNP(chrom, name, position, firstAllele, secondAllele, genoTypes ))
  return(result)
        
with open(args.input, "r") as fin:
  p = re.compile("chr(.+)_\d+_\d+$")
  allSnps = [];
  for line in fin:
    fileInfo = line.split('\t')
    r = p.search(fileInfo[1])
    chrom = int(r.group(1))
    snps = readImpute2File(fileInfo[0], chrom)
    for snp in snps:
      allSnps.append(snp)
      
allSnps = sorted(allSnps, key = lambda x: (x.chrom, x.position))
with open(args.output + ".gen", "w") as fout:
  for snp in allSnps:
    fout.write("%d %s %d %s %s %s\n" % (snp.chrom, snp.name, snp.position, snp.firstAllele, snp.secondAllele, " ".join(snp.genoTypes)))
    
copyfile(args.sample, args.output + ".sample")
os.system("plink2 --data " + args.output + " --ref-first --make-bed --out " + args.output)

