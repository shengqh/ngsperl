import argparse
import sys
import logging
import os
import csv
import gzip
from Bio import SeqIO, bgzf
from subprocess import call

#https://software.broadinstitute.org/gatk/documentation/article.php?id=6926

DEBUG=False
NotDEBUG=not DEBUG

parser = argparse.ArgumentParser(description="Discard the SNV with '*' (spanning deletion) from left trim of GATK4.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i', '--input', action='store', nargs='?', required=NotDEBUG, help='Input VCF file')
parser.add_argument('-o', '--output', action='store', nargs='?', required=NotDEBUG, help='Output VCF file')

args = parser.parse_args()
if DEBUG:
  args.input = "/scratch/cqs/shengq2/macrae_linton/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_gatk4_hc_gvcf_vqsr/result/linton_exomeseq_2118.indels.snp.recal.pass.leftAligned.vcf.gz"
  args.output = "/scratch/cqs/shengq2/macrae_linton/20180913_linton_exomeseq_2118_human_cutadapt/bwa_refine_gatk4_hc_gvcf_vqsr/result/linton_exomeseq_2118.indels.snp.recal.pass.leftAligned.fixed.vcf"

print(args)

logger = logging.getLogger('fixDeletion')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

# logger.info("reading fasta file " + args.fasta)
# fastaMap = {}
# fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')
# for fasta in fasta_sequences:
#   fastaMap[fasta.id] = str(fasta.seq)
#   
logger.info("fixing left aligned deletion ...")
with open(args.output, "w") as fout:
  if args.input.endswith(".gz"):
    fin = gzip.open(args.input, 'rb')
  else:
    fin = open(args.input, "r")
  try:
    while True:
      line = fin.readline()
      if "#CHROM" in line:
        fout.write(line)
        vcfheaders = line.rstrip().split("\t")
        chr_index = 0
        position_index = vcfheaders.index("POS")
        ref_index = vcfheaders.index("REF")
        alt_index = vcfheaders.index("ALT")
        break
      else:
        fout.write(line)
    
    last_locus = ''
    cur_snp_list = []
    for line in fin:
      snv = line.split('\t')
      
      chr = snv[chr_index]
      curposition = snv[position_index]
      cur_locus = chr + ":" + curposition
      
      if cur_locus != last_locus:
        for snp_line in cur_snp_list:
          fout.write(snp_line)
        cur_snp_list = []
        last_locus = cur_locus
      
      if snv[alt_index] == '*':
        continue

      cur_snp_list.append(line)
      
    for snp_line in cur_snp_list:
      fout.write(snp_line)
  finally:
    fin.close()
    
logger.info("done.")
